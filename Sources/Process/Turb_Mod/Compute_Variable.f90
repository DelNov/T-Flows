!==============================================================================!
  subroutine Turb_Mod_Compute_Variable(turb, Sol, curr_dt, ini, phi)
!------------------------------------------------------------------------------!
!   Discretizes and solves transport equations for different turbulent         !
!   variables.                                                                 !
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  type(Turb_Type),   target :: turb
  type(Solver_Type), target :: Sol
  integer, intent(in)       :: curr_dt
  integer, intent(in)       :: ini
  type(Var_Type)            :: phi
!----------------------------------[Locals]------------------------------------!
  type(Field_Type),  pointer :: Flow
  type(Grid_Type),   pointer :: Grid
  type(Var_Type),    pointer :: u, v, w
  type(Var_Type),    pointer :: vis
  real, contiguous,  pointer :: flux(:)
  type(Matrix_Type), pointer :: A
  real, contiguous,  pointer :: b(:)
  integer                    :: s, c, c1, c2
  real                       :: f_ex, f_im
  real                       :: a0, a12, a21
  real                       :: vis_eff
  real                       :: phi_x_f, phi_y_f, phi_z_f
  real                       :: dt
  real                       :: visc_f
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

  call Cpu_Timer % Start('Compute_Turbulence (without solvers)')

  ! Take aliases
  Flow => turb % pnt_flow
  Grid => Flow % pnt_grid
  vis  => turb % vis
  flux => Flow % v_flux % n
  dt   =  Flow % dt
  call Flow % Alias_Momentum(u, v, w)
  call Sol % Alias_Solver      (A, b)

  ! Initialize matrix and right hand side
  A % val(:) = 0.0
  b      (:) = 0.0

  ! Old values (o) and older than old (oo)
  if(ini .eq. 1) then
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
  call Numerics_Mod_Advection_Term(phi, Flow % density, flux, b)

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

    visc_f =        Grid % fw(s)  * Flow % viscosity(c1)   &
           + (1.0 - Grid % fw(s)) * Flow % viscosity(c2)

    vis_eff = visc_f + (    Grid % fw(s)  * turb % vis_t(c1)   &
                     + (1.0-Grid % fw(s)) * turb % vis_t(c2))  &
                     / phi % sigma

    if(turb % model .eq. SPALART_ALLMARAS .or.               &
       turb % model .eq. DES_SPALART)                        &
      vis_eff = visc_f + (    Grid % fw(s)  * vis % n(c1)    &
                       + (1.0-Grid % fw(s)) * vis % n(c2))   &
                       / phi % sigma

    if(turb % model .eq. HYBRID_LES_RANS) then
      vis_eff = visc_f + (    Grid % fw(s)  * turb % vis_t_eff(c1)   &
                       + (1.0-Grid % fw(s)) * turb % vis_t_eff(c2))  &
                       / phi % sigma
    end if
    phi_x_f = Grid % fw(s) * phi % x(c1) + (1.0-Grid % fw(s)) * phi % x(c2)
    phi_y_f = Grid % fw(s) * phi % y(c1) + (1.0-Grid % fw(s)) * phi % y(c2)
    phi_z_f = Grid % fw(s) * phi % z(c1) + (1.0-Grid % fw(s)) * phi % z(c2)

    if(turb % model .eq. K_EPS_ZETA_F    .or.  &
       turb % model .eq. HYBRID_LES_RANS .or.  &
       turb % model .eq. K_EPS) then
      if(c2 < 0 .and. phi % name .eq. 'KIN') then
        if(Grid % Bnd_Cond_Type(c2) .eq. WALL .or.  &
           Grid % Bnd_Cond_Type(c2) .eq. WALLFL) then
          if(turb % y_plus(c1) > 4) then

            phi_x_f = 0.0
            phi_y_f = 0.0
            phi_z_f = 0.0
            vis_eff = 0.0

          end if
        end if
      end if
    end if

    ! Total (exact) diffusive flux
    f_ex = vis_eff * (  phi_x_f * Grid % sx(s)  &
                      + phi_y_f * Grid % sy(s)  &
                      + phi_z_f * Grid % sz(s) )

    a0 = vis_eff * A % fc(s)

    ! Implicit diffusive flux
    f_im = (  phi_x_f * Grid % dx(s)                      &
            + phi_y_f * Grid % dy(s)                      &
            + phi_z_f * Grid % dz(s) ) * a0

    ! Cross diffusion part
    phi % c(c1) = phi % c(c1) + f_ex - f_im
    if(c2  > 0) then
      phi % c(c2) = phi % c(c2) - f_ex + f_im
    end if

    ! Compute coefficients for the sysytem matrix
    a12 = a0
    a21 = a0

    a12 = a12  - min(flux(s), 0.0) * Flow % density(c1)
    a21 = a21  + max(flux(s), 0.0) * Flow % density(c2)

    ! Fill the system matrix
    if(c2  > 0) then
      A % val(A % pos(1,s)) = A % val(A % pos(1,s)) - a12
      A % val(A % dia(c1))  = A % val(A % dia(c1))  + a12
      A % val(A % pos(2,s)) = A % val(A % pos(2,s)) - a21
      A % val(A % dia(c2))  = A % val(A % dia(c2))  + a21
    else if(c2 < 0) then

      ! Outflow is not included because it was causing problems
      if((Grid % Bnd_Cond_Type(c2) .eq. INFLOW)  .or.   &
         (Grid % Bnd_Cond_Type(c2) .eq. WALL)    .or.   &
         (Grid % Bnd_Cond_Type(c2) .eq. PRESSURE).or.   &
         (Grid % Bnd_Cond_Type(c2) .eq. CONVECT) .or.   &
         (Grid % Bnd_Cond_Type(c2) .eq. WALLFL) ) then
        A % val(A % dia(c1)) = A % val(A % dia(c1)) + a12
        b(c1) = b(c1) + a12 * phi % n(c2)
      end if

    end if  ! if c2 < 0

  end do  ! through faces

  ! Cross diffusion terms are treated explicity
  do c = 1, Grid % n_cells
    b(c) = b(c) + phi % c(c)
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
  if(turb % model .eq. K_EPS) then
    if(phi % name .eq. 'KIN') call Turb_Mod_Src_Kin_K_Eps(turb, Sol)
    if(phi % name .eq. 'EPS') call Turb_Mod_Src_Eps_K_Eps(turb, Sol)
    if(Flow % heat_transfer) then
      if(phi % name .eq. 'T2')  call Turb_Mod_Src_T2(turb, Sol)
    end if
  end if

  if(turb % model .eq. K_EPS_ZETA_F .or.  &
     turb % model .eq. HYBRID_LES_RANS) then
    if(phi % name .eq. 'KIN')  call Turb_Mod_Src_Kin_K_Eps_Zeta_F(turb, Sol)
    if(phi % name .eq. 'EPS')  call Turb_Mod_Src_Eps_K_Eps_Zeta_F(turb, Sol)
    if(phi % name .eq. 'ZETA')  &
      call Turb_Mod_Src_Zeta_K_Eps_Zeta_F(turb, Sol, curr_dt)
    if(Flow % heat_transfer) then
      if(phi % name .eq. 'T2')  call Turb_Mod_Src_T2(turb, Sol)
    end if
  end if

  if(turb % model .eq. SPALART_ALLMARAS .or.  &
     turb % model .eq. DES_SPALART) then
    call Turb_Mod_Src_Vis_Spalart_Almaras(turb, Sol)
  end if

  !---------------------------------!
  !                                 !
  !   Solve the equations for phi   !
  !                                 !
  !---------------------------------!

  ! Under-relax the equations
  call Numerics_Mod_Under_Relax(phi, A, b)

  ! Call linear solver to solve the equations
  call Cpu_Timer % Start('Linear_Solver_For_Turbulence')
  call Sol % Bicg(A,              &
                  phi % n,        &
                  b,              &
                  phi % precond,  &
                  phi % mniter,   &
                  phi % eniter,   &
                  phi % tol,      &
                  phi % res)
  call Cpu_Timer % Stop('Linear_Solver_For_Turbulence')

  ! Avoid negative values for all computed turbulent quantities
  do c = 1, Grid % n_cells
    if( phi % n(c) < 0.0 ) phi % n(c) = phi % o(c)
  end do

  ! Set the lower limit of zeta to 1.8
  if(phi % name .eq. 'ZETA') then
    do c = 1, Grid % n_cells
      phi % n(c) = min(phi % n(c), 1.8)
    end do
  end if

  ! Set the lower limit of epsilon 
  if(phi % name .eq. 'EPS'.and.turb % model == HYBRID_LES_RANS) then
    do c = 1, Grid % n_cells
      phi % n(c) = max(phi % n(c), 1.0e-10)
    end do
  end if

  ! Print info on the screen
  if(turb % model .eq. K_EPS        .or.  &
     turb % model .eq. K_EPS_ZETA_F .or.  &
     turb % model .eq. HYBRID_LES_RANS) then
    if(phi % name .eq. 'KIN')  &
      call Info_Mod_Iter_Fill_At(3, 1, phi % name, phi % eniter, phi % res)
    if(phi % name .eq. 'EPS')  &
      call Info_Mod_Iter_Fill_At(3, 2, phi % name, phi % eniter, phi % res)
    if(phi % name .eq. 'ZETA')  &
      call Info_Mod_Iter_Fill_At(3, 3, phi % name, phi % eniter, phi % res)
    if(Flow % heat_transfer) then
      if(phi % name .eq. 'T2')  &
      call Info_Mod_Iter_Fill_At(3, 5, phi % name, phi % eniter, phi % res)
    end if
  end if

  call Flow % Grad_Variable(phi)

  call Cpu_Timer % Stop('Compute_Turbulence (without solvers)')

  end subroutine
