!==============================================================================!
  subroutine Compute_Scalar(Flow, Turb, Vof, Sol, curr_dt, ini, sc)
!------------------------------------------------------------------------------!
!   Purpose: Solve transport equation for user defined scalar.                 !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use User_Mod
  use Work_Mod
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Arguments]--------------------------------!
  type(Field_Type),    target :: Flow
  type(Turb_Type),     target :: Turb
  type(Vof_Type),      target :: Vof
  type(Solver_Type),   target :: Sol
  integer, intent(in)         :: curr_dt
  integer, intent(in)         :: ini
  integer, intent(in)         :: sc
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: Grid
  type(Var_Type),    pointer :: uu, vv, ww, uv, uw, vw
  type(Matrix_Type), pointer :: A
  real, contiguous,  pointer :: b(:)
  type(Face_Type),   pointer :: v_flux
  type(Var_Type),    pointer :: phi
  integer                    :: c, s, c1, c2, row, col
  real                       :: a12, a21
  real                       :: ns, dt
  real                       :: dif_eff, f_ex, f_im
  real                       :: phi_stress, q_exp
  real                       :: phix_f, phiy_f, phiz_f
  real, contiguous,  pointer :: q_turb(:), cross(:)
!------------------------------------------------------------------------------!
!
!  The form of equations which are solved:
!
!     /                /                /
!    |     d phi      |                |
!    | rho ----- dV   | rho u phi dS = | gamma DIV phi dS
!    |      dt        |                |
!   /                /                /
!
!==============================================================================!

  call Profiler % Start('Compute_Scalars (without solvers)')

  call Work % Connect_Real_Cell(q_turb, cross)

  ! Take aliases
  Grid   => Flow % pnt_grid
  v_flux => Flow % v_flux
  phi    => Flow % scalar(sc)
  dt     =  Flow % dt
  call Turb % Alias_Stresses(uu, vv, ww, uv, uw, vw)
  call Sol % Alias_Native   (A, b)

  ! User function
  call User_Mod_Beginning_Of_Compute_Scalar(Flow, Turb, Vof, Sol,  &
                                            curr_dt, ini, sc)

  ! Initialize cross diffusion sources, matrix and right hand side
  cross(:) = 0.0
  A % val(:) = 0.0
  b      (:) = 0.0

  ! Initialize turbulent scalar fluxes
  ! and fluxes coming from interfaces
  q_turb(:) = 0.0

  !--------------------------!
  !   Initialize variables   !
  !--------------------------!

  ! Old values (o and oo)
  if(ini.lt.2) then
    do c = 1, Grid % n_cells
      phi % oo(c) = phi % o(c)
      phi % o (c) = phi % n(c)
    end do
  end if

  ! Gradients
  call Flow % Grad_Variable(phi)

  !---------------!
  !               !
  !   Advection   !
  !               !
  !---------------!
  call Numerics_Mod_Advection_Term(phi, Flow % density, v_flux % n, b)

  !--------------!
  !              !
  !   Difusion   !
  !              !
  !--------------!

  !----------------------------!
  !   Spatial discretization   !
  !----------------------------!
  call Control_Mod_Turbulent_Schmidt_Number(sc_t)  ! get default sc_t (0.9)

  do s = 1, Grid % n_faces

    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    call Turb % Face_Diff_And_Stress(dif_eff, phi_stress, s, sc)

    ! Gradients on the cell face
    phix_f = Grid % fw(s)*phi % x(c1) + (1.0-Grid % fw(s))*phi % x(c2)
    phiy_f = Grid % fw(s)*phi % y(c1) + (1.0-Grid % fw(s))*phi % y(c2)
    phiz_f = Grid % fw(s)*phi % z(c1) + (1.0-Grid % fw(s))*phi % z(c2)

    ! Total (exact) diffusive flux
    f_ex = dif_eff  * (  phix_f * Grid % sx(s)  &
                       + phiy_f * Grid % sy(s)  &
                       + phiz_f * Grid % sz(s))

    ! Implicit diffusive flux
    f_im = dif_eff  * A % fc(s)          &
         * (  phix_f * Grid % dx(s)      &
            + phiy_f * Grid % dy(s)      &
            + phiz_f * Grid % dz(s) )

    ! Calculate the coefficients for the sysytem matrix
    a12 = dif_eff * A % fc(s)
    a21 = dif_eff * A % fc(s)

    a12 = a12  - min(v_flux % n(s), 0.0) * Flow % density(c1)
    a21 = a21  + max(v_flux % n(s), 0.0) * Flow % density(c2)

    ! Cross diffusion part
    cross(c1) = cross(c1) + f_ex - f_im
    if(c2 .gt. 0) then
      cross(c2) = cross(c2) - f_ex + f_im
    end if

    ! Put the influence of turbulent scalar fluxes explicitly in the system
    q_turb(c1) = q_turb(c1) + phi_stress
    if(c2 > 0) then
      q_turb(c2) = q_turb(c2) - phi_stress
    end if

    ! Fill the system matrix
    if(c2 > 0) then
      A % val(A % dia(c1))  = A % val(A % dia(c1)) + a12
      A % val(A % dia(c2))  = A % val(A % dia(c2)) + a21
      A % val(A % pos(1,s)) = A % val(A % pos(1,s)) - a12
      A % val(A % pos(2,s)) = A % val(A % pos(2,s)) - a21
    else if(c2 < 0) then

      ! Outflow is included because of the flux
      ! corrections which also affects velocities
      if( (Var_Mod_Bnd_Cond_Type(phi,c2) .eq. INFLOW) .or.  &
          (Var_Mod_Bnd_Cond_Type(phi,c2) .eq. WALL)   .or.  &
          (Var_Mod_Bnd_Cond_Type(phi,c2) .eq. CONVECT) ) then
        A % val(A % dia(c1)) = A % val(A % dia(c1)) + a12
        b(c1)  = b(c1)  + a12 * phi % n(c2)

      ! In case of wallflux
      else if(Var_Mod_Bnd_Cond_Type(phi,c2) .eq. WALLFL) then
        b(c1) = b(c1) + Grid % s(s) * phi % q(c2)
      end if

    end if

  end do  ! through sides

  !----------------------------------------------!
  !   Explicitly treated diffusion scalar fluxes !
  !   and cross diffusion                        !
  !----------------------------------------------!
  do c = 1, Grid % n_cells

    ! Total explicit heat flux
    q_exp = cross(c) + q_turb(c)

    if(q_exp >= 0) then
      b(c)  = b(c) + q_exp
    else
      A % val(A % dia(c)) = A % val(A % dia(c)) - q_exp / (phi % n(c) + MICRO)
    end if
  end do

  !--------------------!
  !                    !
  !   Inertial terms   !
  !                    !
  !--------------------!

  call Numerics_Mod_Inertial_Term(phi, Flow % density, A, b, dt)

  !-------------------------------------!
  !                                     !
  !   Source terms and wall function    !
  !                                     !
  !-------------------------------------!

  call User_Mod_Source(Flow, phi, A, b)

  !---------------------------------!
  !                                 !
  !   Solve the equations for phi   !
  !                                 !
  !---------------------------------!

  ! Under-relax the equations
  call Numerics_Mod_Under_Relax(phi, A, b)

  call Profiler % Start('Linear_Solver_For_Scalars')

  ! Call linear solver to solve them
  call Sol % Run(phi % solver,     &
                 phi % prec,       &
                 phi % prec_opts,  &
                 A,                &
                 phi % n,          &
                 b,                &
                 phi % mniter,     &
                 phi % eniter,     &
                 phi % tol,        &
                 phi % res)

  call Profiler % Stop('Linear_Solver_For_Scalars')

  read(phi % name(3:4), *) ns  ! reterive the number of scalar
  row = ceiling(ns/6)          ! will be 1 (scal. 1-6), 2 (scal. 6-12), etc.
  col = nint(ns) - (row-1)*6   ! will be in range 1 - 6

  call Info_Mod_Iter_Fill_User_At(row, col, phi % name, phi % eniter, phi % res)

  call Flow % Grad_Variable(phi)

  ! User function
  call User_Mod_End_Of_Compute_Scalar(Flow, Turb, Vof, Sol, curr_dt, ini, sc)

  call Work % Disconnect_Real_Cell(q_turb, cross)

  call Profiler % Stop('Compute_Scalars (without solvers)')

  end subroutine
