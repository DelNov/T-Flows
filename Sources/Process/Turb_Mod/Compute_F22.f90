!==============================================================================!
  subroutine Turb_Mod_Compute_F22(turb, Sol, curr_dt, ini, phi)
!------------------------------------------------------------------------------!
!   Discretizes and solves eliptic relaxation equations for f22.               !
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  type(Turb_Type)           :: turb
  type(Solver_Type), target :: Sol
  integer, intent(in)       :: curr_dt
  integer, intent(in)       :: ini
  type(Var_Type)            :: phi
!----------------------------------[Locals]------------------------------------!
  type(Field_Type),  pointer :: Flow
  type(Grid_Type),   pointer :: Grid
  type(Matrix_Type), pointer :: A
  real, contiguous,  pointer :: b(:)
  integer                    :: s, c, c1, c2
  real                       :: f_ex, f_im
  real                       :: a0, a12, a21
  real                       :: phi_x_f, phi_y_f, phi_z_f
!==============================================================================!
!                                                                              !
!   The form of equations which are solved:                                    !
!                                                                              !
!        /             /              /                                        !
!       |   df22      |   f22        |  f22hg                                  !
!       | - ---- dS + | ------ dV  = | ------- dV                              !
!       |    dy       |  Lsc^2       |  Lsc^2                                  !
!      /             /              /                                          !
!                                                                              !
!   Dimension of the system under consideration                                !
!                                                                              !
!     [A]{f22} = {b}   [kg K/s]                                                !
!                                                                              !
!   Dimensions of certain variables:                                           !
!                                                                              !
!     f22            [1/s]                                                     !
!     Lsc            [m]                                                       !
!------------------------------------------------------------------------------!

  ! Take aliases
  Flow => turb % pnt_flow
  Grid => Flow % pnt_grid
  call Sol % Alias_Native(A, b)

  ! Initialize matrix and right hand side
  A % val(:) = 0.0
  b      (:) = 0.0

  ! Old values (o) and older than old (oo)
  if(ini .eq. 1) then
    do c = 1, Grid % n_cells
      phi % oo(c)   = phi % o(c)
      phi % o (c)   = phi % n(c)
    end do
  end if

  ! New values
  do c = 1, Grid % n_cells
    phi % c(c) = 0.0
  end do

  ! Gradients
  call Flow % Grad_Variable(phi)

  !------------------!
  !                  !
  !     Difusion     !
  !                  !
  !------------------!

  !----------------------------!
  !   Spatial discretization   !
  !----------------------------!
  do s = 1, Grid % n_faces

    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    phi_x_f = Grid % fw(s) * phi % x(c1) + (1.0-Grid % fw(s)) * phi % x(c2)
    phi_y_f = Grid % fw(s) * phi % y(c1) + (1.0-Grid % fw(s)) * phi % y(c2)
    phi_z_f = Grid % fw(s) * phi % z(c1) + (1.0-Grid % fw(s)) * phi % z(c2)

    ! Total (exact) diffusive flux
    f_ex = (  phi_x_f * Grid % sx(s)   &
            + phi_y_f * Grid % sy(s)   &
            + phi_z_f * Grid % sz(s) )

    a0 = A % fc(s)

    ! Implicit diffusive flux
    f_im=(   phi_x_f * Grid % dx(s)        &
           + phi_y_f * Grid % dy(s)        &
           + phi_z_f * Grid % dz(s)) * a0

    ! Cross diffusion part
    phi % c(c1) = phi % c(c1) + f_ex - f_im
    if(c2  > 0) then
      phi % c(c2) = phi % c(c2) - f_ex + f_im
    end if

    ! Calculate the coefficients for the sysytem matrix
    a12 = a0
    a21 = a0

    ! Fill the system matrix
    if(c2  > 0) then
      A % val(A % pos(1,s)) = A % val(A % pos(1,s)) - a12
      A % val(A % dia(c1))  = A % val(A % dia(c1))  + a12
      A % val(A % pos(2,s)) = A % val(A % pos(2,s)) - a21
      A % val(A % dia(c2))  = A % val(A % dia(c2))  + a21
    else if(c2  < 0) then

      ! Inflow
      if( (Grid % Bnd_Cond_Type(c2) .eq. INFLOW)) then
        A % val(A % dia(c1)) = A % val(A % dia(c1)) + a12
        b(c1) = b(c1) + a12 * phi % n(c2)
      end if

      ! Wall and wall flux; solid walls in any case
      if( (Grid % Bnd_Cond_Type(c2) .eq. WALL).or.       &
          (Grid % Bnd_Cond_Type(c2) .eq. WALLFL) ) then
        A % val(A % dia(c1)) = A % val(A % dia(c1)) + a12
        !---------------------------------------------------------------!
        !   Source coefficient is filled in SourceF22.f90 in order to   !
        !   get updated values of f22 on the wall.  Otherwise f22       !
        !   equation does not converge very well                        !
        !   b(c1) = b(c1) + a12 * phi % n(c2)                           !
        !---------------------------------------------------------------!
      end if
    end if

  end do  ! through faces

  ! Cross diffusion terms are treated explicity
  do c = 1, Grid % n_cells
    b(c) = b(c) + phi % c(c)
  end do

  !-------------------------------------!
  !                                     !
  !   Source terms and wall function    !
  !                                     !
  !-------------------------------------!
  if(turb % model .eq. RSM_MANCEAU_HANJALIC) then
    call Turb_Mod_Src_F22_Rsm_Manceau_Hanjalic(turb, Sol)
  else
    call Turb_Mod_Src_F22_K_Eps_Zeta_F(turb, Sol)
  end if

  !---------------------------------!
  !                                 !
  !   Solve the equations for phi   !
  !                                 !
  !---------------------------------!

  ! Underrelax the equations
  call Numerics_Mod_Under_Relax(phi, a, b)

  call Cpu_Timer % Start('Linear_Solver_For_Turbulence')

  ! Call linear solver to solve the equations
  call Sol % Run(phi % solver, phi % prec,  &
                 A,                         &
                 phi % n,                   &
                 b,                         &
                 phi % mniter,              &
                 phi % eniter,              &
                 phi % tol,                 &
                 phi % res)

  call Cpu_Timer % Stop('Linear_Solver_For_Turbulence')

  ! Print info on the screen
  if(turb % model .eq. K_EPS_ZETA_F) then
    call Info_Mod_Iter_Fill_At(3, 4, phi % name, phi % eniter, phi % res)
  else if(turb % model .eq. RSM_MANCEAU_HANJALIC) then
    call Info_Mod_Iter_Fill_At(4, 2, phi % name, phi % eniter, phi % res)
  end if

  call Flow % Grad_Variable(phi)

  end subroutine
