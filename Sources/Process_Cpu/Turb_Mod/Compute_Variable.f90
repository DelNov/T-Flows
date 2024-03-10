!==============================================================================!
  subroutine Compute_Variable(Turb, Sol, phi)
!------------------------------------------------------------------------------!
!   Discretizes and solves transport equations for different turbulent         !
!   variables.                                                                 !
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  class(Turb_Type),  target :: Turb
  type(Solver_Type), target :: Sol
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
  real                       :: visc_f, pr_t1, pr_t2, pr_1, pr_2
  real, contiguous,  pointer :: cross(:)
!------------------------------------------------------------------------------!
!                                                                              !
!  The form of equations which are solved:                                     !
!                                                                              !
!     /               /                /                     /                 !
!    |     dphi      |                | mu_eff              |                  !
!    | rho ---- dV + | rho u phi dS = | ------ DIV phi dS + | G dV             !
!    |      dt       |                |  sigma              |                  !
!   /               /                /                     /                   !
!                                                                              !
!==============================================================================!

  call Profiler % Start('Compute_Variable (without solvers)')

  call Work % Connect_Real_Cell(cross)

  ! Take aliases
  Flow => Turb % pnt_flow
  Grid => Flow % pnt_grid
  vis  => Turb % vis
  flux => Flow % v_flux % n
  dt   =  Flow % dt
  call Flow % Alias_Momentum(u, v, w)
  call Sol % Alias_Native   (A, b)

  ! Initialize cross diffusion sources, matrix and right hand side
  cross(:) = 0.0
  A % val(:) = 0.0
  b      (:) = 0.0

  ! Old values (o) and older than old (oo)
  if(Iter % Current() .eq. 1) then
    do c = Cells_In_Domain_And_Buffers()
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

    vis_eff = visc_f + (    Grid % fw(s)  * Turb % vis_t(c1)   &
                     + (1.0-Grid % fw(s)) * Turb % vis_t(c2))  &
                     / phi % sigma

    if(phi % name .eq. 'T2') then
      pr_t1 = Turb % Prandtl_Turb(c1)
      pr_t2 = Turb % Prandtl_Turb(c2)
      pr_1  = Flow % Prandtl_Numb(c1)
      pr_2  = Flow % Prandtl_Numb(c2)

      visc_f =      Grid % fw(s)  * Flow % viscosity(c1) / pr_1   &
           + (1.0 - Grid % fw(s)) * Flow % viscosity(c2) / pr_2

      vis_eff = visc_f + (    Grid % fw(s)  * Turb % vis_t(c1)  / pr_t1  &
                       + (1.0-Grid % fw(s)) * Turb % vis_t(c2)) / pr_t2  &
                       / phi % sigma
    end if

    if(Turb % model .eq. SPALART_ALLMARAS .or.               &
       Turb % model .eq. DES_SPALART)                        &
      vis_eff = visc_f + (    Grid % fw(s)  * vis % n(c1)    &
                       + (1.0-Grid % fw(s)) * vis % n(c2))   &
                       / phi % sigma

    if(Turb % model .eq. HYBRID_LES_RANS) then
      vis_eff = visc_f + (    Grid % fw(s)  * Turb % vis_t_eff(c1)   &
                       + (1.0-Grid % fw(s)) * Turb % vis_t_eff(c2))  &
                       / phi % sigma
    end if
    phi_x_f = Grid % fw(s) * phi % x(c1) + (1.0-Grid % fw(s)) * phi % x(c2)
    phi_y_f = Grid % fw(s) * phi % y(c1) + (1.0-Grid % fw(s)) * phi % y(c2)
    phi_z_f = Grid % fw(s) * phi % z(c1) + (1.0-Grid % fw(s)) * phi % z(c2)

    if(Turb % model .eq. K_EPS_ZETA_F    .or.  &
       Turb % model .eq. HYBRID_LES_RANS .or.  &
       Turb % model .eq. K_EPS) then
      if(c2 < 0 .and. phi % name .eq. 'KIN') then
        if(Grid % Bnd_Cond_Type(c2) .eq. WALL .or.  &
           Grid % Bnd_Cond_Type(c2) .eq. WALLFL) then
          if(Turb % y_plus(c1) > 4) then

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
    cross(c1) = cross(c1) + f_ex - f_im
    if(c2  > 0) then
      cross(c2) = cross(c2) - f_ex + f_im
    end if

    ! Compute coefficients for the sysytem matrix
    a12 = a0
    a21 = a0

    ! Blend system matrix if desired to do so
    if(phi % blend_matrix) then
      a12 = a12  - min(flux(s), 0.0) * Flow % density(c1)
      a21 = a21  + max(flux(s), 0.0) * Flow % density(c2)
    end if

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
  do c = Cells_In_Domain_And_Buffers()
    b(c) = b(c) + cross(c)
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
  if(Turb % model .eq. K_EPS) then
    if(phi % name .eq. 'KIN') call Turb % Src_Kin_K_Eps(Sol)
    if(phi % name .eq. 'EPS') call Turb % Src_Eps_K_Eps(Sol)
    if(Flow % heat_transfer) then
      if(phi % name .eq. 'T2')  call Turb % Src_T2(Sol)
    end if
  end if

  if(Turb % model .eq. K_EPS_ZETA_F .or.  &
     Turb % model .eq. HYBRID_LES_RANS) then
    if(phi % name .eq. 'KIN')  call Turb % Src_Kin_K_Eps_Zeta_F(Sol)
    if(phi % name .eq. 'EPS')  call Turb % Src_Eps_K_Eps_Zeta_F(Sol)
    if(phi % name .eq. 'ZETA')  &
      call Turb % Src_Zeta_K_Eps_Zeta_F(Sol)
    if(Flow % heat_transfer) then
      if(phi % name .eq. 'T2')  call Turb % Src_T2(Sol)
    end if
  end if

  if(Turb % model .eq. SPALART_ALLMARAS .or.  &
     Turb % model .eq. DES_SPALART) then
    call Turb % Src_Vis_Spalart_Almaras(Sol)
  end if

  !---------------------------------!
  !                                 !
  !   Solve the equations for phi   !
  !                                 !
  !---------------------------------!

  ! Under-relax the equations
  call Numerics_Mod_Under_Relax(phi, A, b)

  call Profiler % Start(String % First_Upper(phi % solver)  //  &
                        ' (solver for turbulence)')

  ! Call linear solver to solve the equations
  call Sol % Run(A, phi, b)

  call Profiler % Stop(String % First_Upper(phi % solver)  //  &
                       ' (solver for turbulence)')

  ! Avoid negative values for all computed turbulent quantities
  do c = Cells_In_Domain_And_Buffers()
    if( phi % n(c) < 0.0 ) phi % n(c) = phi % o(c)
  end do

  ! Set the lower limit of zeta to 1.8
  if(phi % name .eq. 'ZETA') then
    do c = Cells_In_Domain_And_Buffers()
      phi % n(c) = min(phi % n(c), 1.8)
    end do
  end if

  ! Set the lower limit of epsilon 
  if(phi % name .eq. 'EPS') then
    do c = Cells_In_Domain_And_Buffers()
      phi % n(c) = max(phi % n(c), 1.0e-10)
    end do
  end if

  ! Print info on the screen
  if(Turb % model .eq. K_EPS        .or.  &
     Turb % model .eq. K_EPS_ZETA_F .or.  &
     Turb % model .eq. HYBRID_LES_RANS) then
    if(phi % name .eq. 'KIN')  &
      call Info % Iter_Fill_At(3, 1, phi % name, phi % res, phi % niter)
    if(phi % name .eq. 'EPS')  &
      call Info % Iter_Fill_At(3, 2, phi % name, phi % res, phi % niter)
    if(phi % name .eq. 'ZETA')  &
      call Info % Iter_Fill_At(3, 3, phi % name, phi % res, phi % niter)
    if(Flow % heat_transfer) then
      if(phi % name .eq. 'T2')  &
      call Info % Iter_Fill_At(3, 5, phi % name, phi % res, phi % niter)
    end if
  end if

  call Flow % Grad_Variable(phi)

  call Work % Disconnect_Real_Cell(cross)

  call Profiler % Stop('Compute_Variable (without solvers)')

  end subroutine
