!==============================================================================!
  subroutine Compute_Energy(Process, Flow, Turb, Vof, Sol)
!------------------------------------------------------------------------------!
!   Purpose: Solve transport equation for scalar (such as temperature)         !
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Arguments]--------------------------------!
  class(Process_Type)         :: Process
  type(Field_Type),    target :: Flow
  type(Turb_Type),     target :: Turb
  type(Vof_Type),      target :: Vof
  type(Solver_Type),   target :: Sol
!-----------------------------------[Locals]-----------------------------------! 
  type(Grid_Type),   pointer :: Grid
  type(Var_Type),    pointer :: u, v, w, t
  type(Face_Type),   pointer :: v_flux
  type(Matrix_Type), pointer :: A
  real, contiguous,  pointer :: b(:)
  integer                    :: c, s, c1, c2
  real                       :: a12, a21, con_eff, dt
  real                       :: f_ex, f_im, tx_f, ty_f, tz_f, t_stress, q_exp
  real, contiguous,  pointer :: cap_dens(:), q_int(:), q_turb(:), cross(:)
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Process)
!------------------------------------------------------------------------------!
!
!  The form of equations which are solved:
!
!     /                /                 /
!    |        dT      |                 |
!    | rho Cp -- dV   | rho u Cp T dS = | lambda  DIV T dS
!    |        dt      |                 |
!   /                /                 /
!
!   [A]{T} = {b}
!
!------------------------------------------------------------------------------!
!   Dimensions of certain variables:                                           !
!                                                                              !
!   lambda <-> conductivity, con_w
!   rho    <-> density
!   Cp     <-> capacity
!   T      <-> t % n, t % o, t % oo
!   heat capacity               capacity          [J/(kg K)]
!   thermal conductivity        conductivity      [W/(m K)] ([W = J/s])
!   density                     density           [kg/m^3]
!   flux                        v_flux            [m^3/s]
!   left  hand s.               A                 [J/(s K)]
!   temperature                 t % n             [K]
!   right hand s.               b                 [J/s]
!   Turb. thermal conductivity  con_t_f           [W/(m K)]
!   Turb. hear flux             ut_x_cap_dens_s   [J/(m^2 s)]
!   Turb. stress                t_stress          [J/s]
!   Turb. viscosity             vis_t             [kg/(m s)]
!
!   User_Mod variables:
!   bulk flux                   bulk % flux_x     [kg/s]
!   heat flux                   Flow % heat_flux  [W/m^2]
!==============================================================================!

  if(.not. Flow % heat_transfer) return

  call Profiler % Start('Compute_Energy (without solvers)')

  call Work % Connect_Real_Cell(cap_dens, q_int, q_turb, cross)

  ! Take aliases
  Grid   => Flow % pnt_grid
  v_flux => Flow % v_flux
  dt     =  Flow % dt
  call Flow % Alias_Momentum(u, v, w)
  call Flow % Alias_Energy  (t)
  call Sol % Alias_Native   (A, b)

  ! User function
  call User_Mod_Beginning_Of_Compute_Energy(Flow, Turb, Vof, Sol)

  ! Initialize cross diffusion sources, matrix and right hand side
  cross  (:) = 0.0
  A % val(:) = 0.0
  b      (:) = 0.0

  ! Initialize turbulent heat fluxes (used with RSM models)
  ! and fluxes coming from interfaces (cases with boiling)
  q_turb(:) = 0.0
  q_int (:) = 0.0

  ! Old values (o and oo)
  if(Iter % Current() .eq. 1) then
    do c = 1, Grid % n_cells
      t % oo(c) = t % o(c)
      t % o (c) = t % n(c)
    end do
  end if

  ! Gradients
  if(Flow % mass_transfer_model .eq. NO_MASS_TRANSFER) then
    call Flow % Grad_Variable(t)

  ! If mass transfer, estimate the mass transfer due to heat fluxes,
  ! which also compute gradients with saturation temperature at interface
  else
    call Vof % Mass_Transfer_Estimate()
  end if

  ! Compute helping variable
  do c = -Grid % n_bnd_cells, Grid % n_cells
    cap_dens(c) = Flow % capacity(c) * Flow % density(c)
  end do

  !---------------!
  !               !
  !   Advection   !
  !               !
  !---------------!
  call Numerics_Mod_Advection_Term(t, cap_dens, v_flux % n, b)

  !--------------!
  !              !
  !   Difusion   !
  !              !
  !--------------!

  !----------------------------!
  !   Spatial discretization   !
  !----------------------------!
  do s = 1, Grid % n_faces

    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    call Turb % Face_Cond_And_Stress(con_eff, t_stress, s)

    if(Flow % mass_transfer_model .ne. NO_MASS_TRANSFER) then
      if(Vof % fun % n(c1) < 0.5 .and.  &
         Vof % fun % n(c2) < 0.5) con_eff = Vof % phase_cond(0)
      if(Vof % fun % n(c1) > 0.5 .and.  &
         Vof % fun % n(c2) > 0.5) con_eff = Vof % phase_cond(1)
    end if

    ! Gradients on the cell face (fw corrects situation close to the wall)
    tx_f = Grid % fw(s) * t % x(c1) + (1.0-Grid % fw(s)) * t % x(c2)
    ty_f = Grid % fw(s) * t % y(c1) + (1.0-Grid % fw(s)) * t % y(c2)
    tz_f = Grid % fw(s) * t % z(c1) + (1.0-Grid % fw(s)) * t % z(c2)

    ! Total (exact) diffusion flux
    ! This is last term in equation 2.33 in Denner's thesis because:
    ! Grid % sx (T-Flows) == n_f * A_f (Denner)
    f_ex = con_eff * (   tx_f * Grid % sx(s)   &
                       + ty_f * Grid % sy(s)   &
                       + tz_f * Grid % sz(s))

    ! Implicit part of the diffusion flux, treated by linear system
    ! This is also term in equation 2.33 in Denner's thesis because:
    ! A % fc * Grid % dx (T-Flows) == alpha_f * s_f * A_f (Denner)
    f_im = con_eff * A % fc(s) * (   tx_f * Grid % dx(s)    &
                                   + ty_f * Grid % dy(s)    &
                                   + tz_f * Grid % dz(s) )


    !------------------------------------------------------!
    !                                                      !
    !   Calculate the coefficients for the system matrix   !
    !                                                      !
    !------------------------------------------------------!
    a12 = con_eff * A % fc(s)
    a21 = con_eff * A % fc(s)

    ! Blend system matrix if desired to do so
    if(t % blend_matrix) then
      a12 = a12 - min(v_flux % n(s), 0.0)  &
                   * Flow % capacity(c1) * Flow % density(c1)  ! Flow: 1 -> 2
      a21 = a21 + max(v_flux % n(s), 0.0)  &
                   * Flow % capacity(c2) * Flow % density(c2)  ! Flow: 2 -> 1
    end if

    !-----------------------------------------------------!
    !   In case of mass transfer, detach the two phases   !
    !      and add heat transferred to the interface      !
    !-----------------------------------------------------!
    if(Flow % mass_transfer_model == TEMPERATURE_GRADIENTS) then
      if(Vof % Front % intersects_face(s)) then
        a12  = 0.0
        a21  = 0.0
        f_ex = 0.0
        f_im = 0.0
        A % val(A % dia(c1)) = A % val(A % dia(c1)) + Vof % a12(s)
        b(c1) = b(c1) + Vof % a12(s) * Vof % t_sat
        A % val(A % dia(c2)) = A % val(A % dia(c2)) + Vof % a21(s)
        b(c2) = b(c2) + Vof % a21(s) * Vof % t_sat
      end if
    end if

    ! Cross diffusion part
    cross(c1) = cross(c1) + f_ex - f_im
    if(c2 > 0) then
      cross(c2) = cross(c2) - f_ex + f_im
    end if

    ! Put the influence of turbulent heat fluxes explicitly in the system
    q_turb(c1) = q_turb(c1) + t_stress
    if(c2 > 0) then
      q_turb(c2) = q_turb(c2) - t_stress
    end if

    !----------------------------!
    !   Fill the system matrix   !
    !----------------------------!
    if(c2 > 0) then
      A % val(A % dia(c1))  = A % val(A % dia(c1)) + a12
      A % val(A % dia(c2))  = A % val(A % dia(c2)) + a21
      A % val(A % pos(1,s)) = A % val(A % pos(1,s)) - a12
      A % val(A % pos(2,s)) = A % val(A % pos(2,s)) - a21
    else if(c2 .lt. 0) then
      ! Outflow is included because of the flux
      ! corrections which also affects velocities
      if( (Var_Mod_Bnd_Cond_Type(t, c2) .eq. INFLOW) .or.  &
          (Var_Mod_Bnd_Cond_Type(t, c2) .eq. WALL)   .or.  &
          (Var_Mod_Bnd_Cond_Type(t, c2) .eq. CONVECT) ) then
        A % val(A % dia(c1)) = A % val(A % dia(c1)) + a12
        b(c1)  = b(c1)  + a12 * t % n(c2)
      ! In case of wallflux 
      else if(Var_Mod_Bnd_Cond_Type(t, c2) .eq. WALLFL) then
        b(c1) = b(c1) + Grid % s(s) * t % q(c2)
      end if
    end if

  end do  ! through sides

  ! heat sink or source due to mass transfer
  if(Flow % mass_transfer_model == LEE) then

    ! Heat sink or source for Lee model (Yohei)
    do c = Cells_In_Domain_And_Buffers()
      ! Units: (kg/s) * (J/kg) = W
      b(c) = b(c) - Vof % m_dot(c) * Vof % latent_heat
    end do
  end if

  !----------------------------------------------!
  !   Explicitly treated diffusion heat fluxes   !
  !   cross diffusion, and heat from interface   !
  !----------------------------------------------!
  do c = 1, Grid % n_cells

    ! Total explicit heat flux
    q_exp = cross(c) + q_turb(c) + q_int(c)

    if(q_exp >= 0) then
      b(c)  = b(c) + q_exp
    else
      A % val(A % dia(c)) = A % val(A % dia(c)) - q_exp / (t % n(c) + MICRO)
    end if
  end do

  !--------------------!
  !                    !
  !   Inertial terms   !
  !                    !
  !--------------------!
  do c = -Grid % n_bnd_cells, Grid % n_cells
    cap_dens(c) = Flow % capacity(c) * Flow % density(c)
    if(Flow % mass_transfer_model .ne. NO_MASS_TRANSFER) then
      if(Vof % fun % n(c) > 0.5) then
        cap_dens(c) = Vof % phase_capa(1) * Vof % phase_dens(1)
      else if(Vof % fun % n(c) < 0.5) then
        cap_dens(c) = Vof % phase_capa(0) * Vof % phase_dens(0)
      end if
    end if
  end do
  call Numerics_Mod_Inertial_Term(t, cap_dens, A, b, dt)

  !--------------------!
  !                    !
  !   User source(s)   !
  !                    !
  !--------------------!
  call User_Mod_Source(Flow, t, A, b)

  !-------------------------------!
  !                               !
  !   Solve the equations for t   !
  !                               !
  !-------------------------------!

  ! Under-relax the equations
  call Numerics_Mod_Under_Relax(t, A, b)

  call Profiler % Start(String % First_Upper(t % solver)  //  &
                        ' (solver for energy)')

  ! Call linear solver to solve the equations
  call Sol % Run(A, t, b)

  call Profiler % Stop(String % First_Upper(t % solver)  //  &
                       ' (solver for energy)')

  ! Print some info on the screen
  call Info % Iter_Fill_At(1, 6, t % name, t % res, t % niter)

  ! Gradients
  if(Flow % mass_transfer_model .eq. NO_MASS_TRANSFER) then
    call Flow % Grad_Variable(t)

  ! If mass transfer, compute gradients with saturation temperature
  else
    call Vof % Grad_Variable_With_Front(t, Vof % t_sat)
  end if

  ! User function
  call User_Mod_End_Of_Compute_Energy(Flow, Turb, Vof, Sol)

  call Work % Disconnect_Real_Cell(cap_dens, q_int, q_turb, cross)

  call Profiler % Stop('Compute_Energy (without solvers)')

  end subroutine
