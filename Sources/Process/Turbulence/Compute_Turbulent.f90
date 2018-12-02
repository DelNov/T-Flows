!==============================================================================!
  subroutine Compute_Turbulent(flow, sol, dt, ini, phi, n_step)
!------------------------------------------------------------------------------!
!   Discretizes and solves transport equations for different turbulent         !
!   variables.                                                                 !
!------------------------------------------------------------------------------!
!---------------------------------[Modules]------------------------------------!
  use Field_Mod
  use Les_Mod
  use Rans_Mod
  use Comm_Mod
  use Var_Mod,      only: Var_Type
  use Grid_Mod,     only: Grid_Type
  use Grad_Mod
  use Info_Mod,     only: Info_Mod_Iter_Fill_At
  use Control_Mod
  use Numerics_Mod, only: CENTRAL, LINEAR, PARABOLIC
  use Solver_Mod,   only: Solver_Type, Bicg, Cg, Cgs
  use Matrix_Mod,   only: Matrix_Type
  use Work_Mod,     only: phi_x   => r_cell_01,  &
                          phi_y   => r_cell_02,  &
                          phi_z   => r_cell_03,  &
                          phi_min => r_cell_04,  &
                          phi_max => r_cell_05
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  type(Field_Type),  target :: flow
  type(Solver_Type), target :: sol
  real                      :: dt
  integer                   :: ini
  type(Var_Type)            :: phi
  integer                   :: n_step
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type),   pointer :: grid
  type(Var_Type),    pointer :: u, v, w
  real,              pointer :: flux(:)
  type(Matrix_Type), pointer :: a
  real,              pointer :: b(:)
  integer                    :: s, c, c1, c2, niter
  real                       :: f_ex, f_im
  real                       :: phis
  real                       :: a0, a12, a21
  real                       :: ini_res, tol
  real                       :: vis_eff
  real                       :: phi_x_f, phi_y_f, phi_z_f
  character(len=80)          :: precond
  integer                    :: adv_scheme  ! advection scheme
  real                       :: blend       ! blending (1.0 central; 0.0 upwind)
  integer                    :: td_scheme   ! time-disretization for inerita
  real                       :: urf         ! under-relaxation factor
!==============================================================================!
!                                                                              !
!  The form of equations which are solved:                                     !
!                                                                              !
!     /               /                /                     /                 !
!    |     dphi      |                | mu_eff              |                  !
!    | rho ---- dV + | rho u phi dS = | ------ DIV phi dS + | G dV             !
!    |      dt       |                |  sigma              |                  !
!   /               /                /                     /                   !
!                                                                              !
!------------------------------------------------------------------------------!

  ! Take aliases
  grid => flow % pnt_grid
  flux => flow % flux
  u    => flow % u
  v    => flow % v
  w    => flow % w
  a    => sol  % a
  b    => sol  % b % val

  ! Initialize matrix and right hand side
  a % val(:) = 0.0
  b      (:) = 0.0

  ! Old values (o) and older than old (oo)
  if(ini .eq. 1) then
    do c = 1, grid % n_cells
      phi % oo(c) = phi % o(c)
      phi % o (c) = phi % n(c)
    end do
  end if

  ! Gradients
  call Grad_Mod_For_Phi(grid, phi % n, 1, phi_x, .true.)
  call Grad_Mod_For_Phi(grid, phi % n, 2, phi_y, .true.)
  call Grad_Mod_For_Phi(grid, phi % n, 3, phi_z, .true.)

  !---------------!
  !               !
  !   Advection   !
  !               !
  !---------------!

  ! Retreive advection scheme and blending coefficient
  call Control_Mod_Advection_Scheme_For_Turbulence(adv_scheme)
  call Control_Mod_Blending_Coefficient_For_Turbulence(blend)

  ! Compute phimax and phimin
  if(adv_scheme .ne. CENTRAL) then
    call Calculate_Minimum_Maximum(grid, phi % n, phi_min, phi_max)
    goto 1
  end if

  ! New values
1 do c = 1, grid % n_cells
    phi % a(c) = 0.0
    phi % c(c) = 0.0
  end do

  !----------------------------!
  !   Spatial Discretization   !
  !----------------------------!
  do s = 1, grid % n_faces

    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    ! Velocities on "orthogonal" cell centers
    if(c2 > 0) then
      phis =        grid % f(s)  * phi % n(c1)  &
           + (1.0 - grid % f(s)) * phi % n(c2)

      ! Compute phis with desired advection scheme
      if(adv_scheme .ne. CENTRAL) then
        call Advection_Scheme(flow, phis, s, phi % n, phi_min, phi_max,  &
                              phi_x, phi_y, phi_z,                       &
                              grid % dx, grid % dy, grid % dz,           &
                              adv_scheme, blend)
      end if

      ! Compute advection term
      if(c2  > 0) then
        phi % a(c1) = phi % a(c1)-flux(s) * phis
        phi % a(c2) = phi % a(c2)+flux(s) * phis
      else
        phi % a(c1) = phi % a(c1)-flux(s) * phis
      end if

      ! Store upwinded part of the advection term in "c"
      if(flux(s)  < 0) then   ! from c2 to c1
      phi % c(c1) = phi % c(c1) - flux(s) * phi % n(c2)
        if(c2  > 0) then
          phi % c(c2) = phi % c(c2) + flux(s) * phi % n(c2)
        end if
      else
        phi % c(c1) = phi % c(c1) - flux(s) * phi % n(c1)
        if(c2  > 0) then
          phi % c(c2) = phi % c(c2) + flux(s) * phi % n(c1)
        end if
      end if
    end if     ! c2 > 0
  end do    ! through sides

  !------------------------------------------------!
  !   Source term contains difference between      !
  !   explicity and implicitly treated advection   !
  !------------------------------------------------!
  do c = 1, grid % n_cells
    b(c) = b(c) + (phi % a(c) - phi % c(c))
  end do

  !------------------!
  !                  !
  !     Difusion     !
  !                  !
  !------------------!

  ! Set c terms back to zero
  do c = 1, grid % n_cells
    phi % c(c) = 0.0
  end do

  !----------------------------!
  !   Spatial discretization   !
  !----------------------------!
  do s = 1, grid % n_faces

    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    vis_eff = viscosity + (    grid % fw(s)  * vis_t(c1)         &
                        + (1.0-grid % fw(s)) * vis_t(c2))        &
                        / phi % sigma

    if(turbulence_model .eq. SPALART_ALLMARAS .or.               &
       turbulence_model .eq. DES_SPALART)                        &
      vis_eff = viscosity + (    grid % fw(s)  * vis % n(c1)     &
                          + (1.0-grid % fw(s)) * vis % n(c2))    &
                          / phi % sigma

    if(turbulence_model .eq. HYBRID_LES_RANS) then
      vis_eff = viscosity + (    grid % fw(s)  * vis_t_eff(c1)   &
                          + (1.0-grid % fw(s)) * vis_t_eff(c2))  &
                          / phi % sigma
    end if
    phi_x_f = grid % fw(s) * phi_x(c1) + (1.0-grid % fw(s)) * phi_x(c2)
    phi_y_f = grid % fw(s) * phi_y(c1) + (1.0-grid % fw(s)) * phi_y(c2)
    phi_z_f = grid % fw(s) * phi_z(c1) + (1.0-grid % fw(s)) * phi_z(c2)

    if(turbulence_model .eq. K_EPS_ZETA_F    .or.  &
       turbulence_model .eq. HYBRID_LES_RANS .or.  &
       turbulence_model .eq. K_EPS) then
      if(c2 < 0 .and. phi % name .eq. 'KIN') then
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL .or.  &
           Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then
          if(y_plus(c1) > 4) then

            phi_x_f = 0.0
            phi_y_f = 0.0
            phi_z_f = 0.0
            vis_eff = 0.0

          end if
        end if
      end if
    end if

    ! Total (exact) diffusive flux
    f_ex = vis_eff * (  phi_x_f * grid % sx(s)  &
                     + phi_y_f * grid % sy(s)  &
                     + phi_z_f * grid % sz(s) )

    a0 = vis_eff * a % fc(s)

    ! Implicit diffusive flux
    ! (this is a very crude approximation: f_coef is
    !  not corrected at interface between materials)
    f_im = (  phi_x_f * grid % dx(s)                      &
            + phi_y_f * grid % dy(s)                      &
            + phi_z_f * grid % dz(s) ) * a0

    ! Cross diffusion part
    phi % c(c1) = phi % c(c1) + f_ex - f_im
    if(c2  > 0) then
      phi % c(c2) = phi % c(c2) - f_ex + f_im
    end if

    ! Compute coefficients for the sysytem matrix
    a12 = a0
    a21 = a0

    a12 = a12  - min(flux(s), real(0.0))
    a21 = a21  + max(flux(s), real(0.0))

    ! Fill the system matrix
    if(c2  > 0) then
      a % val(a % pos(1,s)) = a % val(a % pos(1,s)) - a12
      a % val(a % dia(c1))  = a % val(a % dia(c1))  + a12
      a % val(a % pos(2,s)) = a % val(a % pos(2,s)) - a21
      a % val(a % dia(c2))  = a % val(a % dia(c2))  + a21
    else if(c2  < 0) then

      ! Outflow is not included because it was causing problems
      if((Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. INFLOW)  .or.                 &
         (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL)    .or.                 &
         (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. PRESSURE).or.                 &
         (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. CONVECT) .or.                 &
         (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) ) then
        a % val(a % dia(c1)) = a % val(a % dia(c1)) + a12
        b(c1) = b(c1) + a12 * phi % n(c2)
      end if
    end if

  end do  ! through faces

  ! Cross diffusion terms are treated explicity
  do c = 1, grid % n_cells
    b(c) = b(c) + phi % c(c)
  end do

  !--------------------!
  !                    !
  !   Inertial terms   !
  !                    !
  !--------------------!

  call Control_Mod_Time_Integration_Scheme(td_scheme)

  ! Two time levels; linear interpolation
  if(td_scheme .eq. LINEAR) then
    do c = 1, grid % n_cells
      a0 = density*grid % vol(c)/dt
      a % val(a % dia(c)) = a % val(a % dia(c)) + a0
      b(c) = b(c) + a0 * phi % o(c)
    end do
  end if

  ! Three time levels; parabolic interpolation
  if(td_scheme .eq. PARABOLIC) then
    do c = 1, grid % n_cells
      a0 = density*grid % vol(c)/dt
      a % val(a % dia(c)) = a % val(a % dia(c)) + 1.5 * a0
      b(c) = b(c) + 2.0*a0 * phi % o(c) - 0.5*a0 * phi % oo(c)
    end do
  end if

  !-------------------------------------!
  !                                     !
  !   Source terms and wall function    !
  !   (Check if it is good to call it   !
  !    before the under relaxation ?)   !
  !                                     !
  !-------------------------------------!
  if(turbulence_model .eq. K_EPS) then
    if(phi % name .eq. 'KIN') call Source_Kin_K_Eps(flow, sol)
    if(phi % name .eq. 'EPS') call Source_Eps_K_Eps(flow, sol)
  end if

  if(turbulence_model .eq. K_EPS_ZETA_F .or.  &
     turbulence_model .eq. HYBRID_LES_RANS) then
    if(phi % name .eq. 'KIN')  call Source_Kin_K_Eps_Zeta_F(flow, sol)
    if(phi % name .eq. 'EPS')  call Source_Eps_K_Eps_Zeta_F(flow, sol)
    if(phi % name .eq. 'ZETA') call Source_Zeta_K_Eps_Zeta_F(grid, sol, n_step)
  end if

  if(turbulence_model .eq. SPALART_ALLMARAS .or.  &
     turbulence_model .eq. DES_SPALART) then
    call Source_Vis_Spalart_Almaras(grid, sol, phi_x, phi_y, phi_z)
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
    b(c) = b(c) + a % val(a % dia(c)) * (1.0 - urf) * phi % n(c) / urf
    a % val(a % dia(c)) = a % val(a % dia(c)) / urf
  end do

  ! Get tolerance for linear solvers
  call Control_Mod_Tolerance_For_Turbulence_Solver(tol)

  ! Get matrix precondioner
  call Control_Mod_Preconditioner_For_System_Matrix(precond)

  ! Set the number of iterations then overwrite with control file if specified
  niter = 6
  call Control_Mod_Max_Iterations_For_Turbulence_Solver(niter)

  call Bicg(sol,      &
            phi % n,  &
            b,        &
            precond,  &
            niter,    &
            tol,      &
            ini_res,  &
            phi % res)

  do c = 1, grid % n_cells
    if( phi % n(c) < 0.0 ) phi % n(c) = phi % o(c)
  end do

  if(turbulence_model .eq. K_EPS        .or.  &
     turbulence_model .eq. K_EPS_ZETA_F .or.  &
     turbulence_model .eq. HYBRID_LES_RANS) then
    if(phi % name .eq. 'KIN')  &
      call Info_Mod_Iter_Fill_At(3, 1, phi % name, niter, phi % res)
    if(phi % name .eq. 'EPS')  &
      call Info_Mod_Iter_Fill_At(3, 2, phi % name, niter, phi % res)
    if(phi % name .eq. 'ZETA')  &
      call Info_Mod_Iter_Fill_At(3, 3, phi % name, niter, phi % res)
  end if

  call Comm_Mod_Exchange_Real(grid, phi % n)

  end subroutine
