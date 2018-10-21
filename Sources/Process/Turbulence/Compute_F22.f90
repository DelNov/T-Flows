!==============================================================================!
  subroutine Compute_F22(grid, ini, phi)
!------------------------------------------------------------------------------!
!   Discretizes and solves eliptic relaxation equations for f22.               !
!------------------------------------------------------------------------------!
!---------------------------------[Modules]------------------------------------!
  use Flow_Mod
  use Les_Mod
  use Rans_Mod
  use Comm_Mod
  use Var_Mod
  use Grid_Mod
  use Grad_Mod
  use Info_Mod
  use Numerics_Mod
  use Solvers_Mod, only: Bicg, Cg, Cgs
  use Control_Mod
  use Work_Mod,    only: phi_x => r_cell_01,  &
                         phi_y => r_cell_02,  &
                         phi_z => r_cell_03           
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  type(Grid_Type) :: grid
  integer         :: ini
  type(Var_Type)  :: phi
!----------------------------------[Locals]------------------------------------!
  integer           :: s, c, c1, c2, niter
  real              :: f_ex, f_im 
  real              :: a0, a12, a21
  real              :: ini_res, tol
  real              :: phi_x_f, phi_y_f, phi_z_f
  character(len=80) :: precond
  real              :: urf           ! under-relaxation factor
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

  a % val = 0.0

  b = 0.0

  ! Old values (o) and older than old (oo)
  if(ini .eq. 1) then
    do c = 1, grid % n_cells
      phi % oo(c)   = phi % o(c)
      phi % o (c)   = phi % n(c)
    end do
  end if

  ! New values
  do c = 1, grid % n_cells
    phi % c(c) = 0.0
  end do

  ! Gradients
  call Grad_Mod_For_Phi(grid, phi % n, 1, phi_x, .true.)
  call Grad_Mod_For_Phi(grid, phi % n, 2, phi_y, .true.)
  call Grad_Mod_For_Phi(grid, phi % n, 3, phi_z, .true.)

  !------------------!
  !                  !
  !     Difusion     !
  !                  !
  !------------------!

  !----------------------------!
  !   Spatial discretization   !
  !----------------------------!
  do s = 1, grid % n_faces

    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    phi_x_f = fw(s) * phi_x(c1) + (1.0 - fw(s)) * phi_x(c2)
    phi_y_f = fw(s) * phi_y(c1) + (1.0 - fw(s)) * phi_y(c2)
    phi_z_f = fw(s) * phi_z(c1) + (1.0 - fw(s)) * phi_z(c2)


    ! Total (exact) diffusive flux
    f_ex = (  phi_x_f * grid % sx(s)   &
            + phi_y_f * grid % sy(s)   &
            + phi_z_f * grid % sz(s) )

    a0 = f_coef(s)

    ! Implicit diffusive flux
    f_im=(   phi_x_f * grid % dx(s)        &
           + phi_y_f * grid % dy(s)        &
           + phi_z_f * grid % dz(s)) * a0

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
      a % val(a % pos(1,s)) = a % val(a % pos(1,s)) - a12
      a % val(a % dia(c1))  = a % val(a % dia(c1))  + a12
      a % val(a % pos(2,s)) = a % val(a % pos(2,s)) - a21
      a % val(a % dia(c2))  = a % val(a % dia(c2))  + a21
    else if(c2  < 0) then

      ! Inflow
      if( (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. INFLOW)) then
        a % val(a % dia(c1)) = a % val(a % dia(c1)) + a12
        b(c1) = b(c1) + a12 * phi % n(c2)
      end if

      ! Wall and wall flux; solid walls in any case
      if( (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL).or.       &
          (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) ) then
        a % val(a % dia(c1)) = a % val(a % dia(c1)) + a12
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
  do c = 1, grid % n_cells
    b(c) = b(c) + phi % c(c)
  end do

  !-------------------------------------!
  !                                     !
  !   Source terms and wall function    !
  !   (Check if it is good to call it   !
  !    before the under relaxation ?)   !
  !                                     !
  !-------------------------------------!
  if(turbulence_model .eq. RSM_MANCEAU_HANJALIC) then
    call Source_F22_Rsm_Manceau_Hanjalic(grid)
  else
    call Source_F22_K_Eps_Zeta_F(grid)
  end if

  !---------------------------------!
  !                                 !
  !   Solve the equations for phi   !
  !                                 !
  !---------------------------------!

  ! Set under-relaxation factor then overwrite with control file if specified
  urf = 1.0
  call Control_Mod_Simple_Underrelaxation_For_Turbulence(urf)

  do c = 1, grid % n_cells
    b(c) = b(c) + a % val(a % dia(c)) * (1.0 - urf)*phi % n(c) / urf
    a % val(a % dia(c)) = a % val(a % dia(c)) / urf
  end do 

  ! Get tolerance for linear solver
  call Control_Mod_Tolerance_For_Turbulence_Solver(tol)

  ! Get matrix precondioner
  call Control_Mod_Preconditioner_For_System_Matrix(precond)

  ! Set the number of iterations then overwrite with control file if specified
  niter = 6
  call Control_Mod_Max_Iterations_For_Turbulence_Solver(niter)

  call Cg(a, phi % n, b, precond, niter, tol, ini_res, phi % res)

  call Info_Mod_Iter_Fill_At(3, 4, phi % name, niter, phi % res)

  call Comm_Mod_Exchange_Real(grid, phi % n)

  end subroutine
