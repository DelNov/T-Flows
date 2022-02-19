!==============================================================================!
  subroutine Compute_Energy(Flow, turb, Vof, Sol, curr_dt, ini)
!------------------------------------------------------------------------------!
!   Purpose: Solve transport equation for scalar (such as temperature)         !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use User_Mod
  use Work_Mod, only: cap_dens => r_cell_01,  &
                      q_int    => r_cell_02,  &
                      q_turb   => r_cell_03      ! turbulent heat fluxes
!------------------------------------------------------------------------------!
!   When using Work_Mod, calling sequence should be outlined                   !
!                                                                              !
!   Main_Pro                (allocates Work_Mod)                               !
!     |                                                                        |
!     +----> Compute_Energy (safe to use r_cell_01..03)                        !
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Arguments]--------------------------------!
  type(Field_Type),    target :: Flow
  type(Turb_Type),     target :: turb
  type(Vof_Type),      target :: Vof
  type(Solver_Type),   target :: Sol
  integer, intent(in)         :: curr_dt
  integer, intent(in)         :: ini
!-----------------------------------[Locals]-----------------------------------! 
  type(Grid_Type),   pointer :: Grid
  type(Var_Type),    pointer :: u, v, w, t
  type(Var_Type),    pointer :: ut, vt, wt
  type(Face_Type),   pointer :: v_flux
  type(Matrix_Type), pointer :: A
  real, contiguous,  pointer :: b(:)
  integer                    :: c, s, c1, c2
  real                       :: a12, a21, con_eff, dt
  real                       :: f_ex, f_im, tx_f, ty_f, tz_f, t_stress, q_exp
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
!   turb. thermal conductivity  con_t_f           [W/(m K)]
!   turb. hear flux             ut                [(m K)/s]
!   turb. hear flux             ut_x_cap_dens_s   [J/(m^2 s)]
!   turb. stress                t_stress          [J/s]
!   turb. viscosity             vis_t             [kg/(m s)]
!
!   User_Mod variables:
!   bulk flux                   bulk % flux_x     [kg/s]
!   heat flux                   Flow % heat_flux  [W/m^2]
!==============================================================================!

  call Cpu_Timer % Start('Compute_Energy (without solvers)')

  ! Take aliases
  Grid   => Flow % pnt_grid
  v_flux => Flow % v_flux
  dt     =  Flow % dt
  call Flow % Alias_Momentum(u, v, w)
  call Flow % Alias_Energy  (t)
  call Sol % Alias_Native   (A, b)

  ! User function
  call User_Mod_Beginning_Of_Compute_Energy(Flow, turb, Vof, Sol, curr_dt, ini)

  ! Initialize advection and cross diffusion sources, matrix and right hand side
  t % a  (:) = 0.0
  t % c  (:) = 0.0
  A % val(:) = 0.0
  b      (:) = 0.0

  ! Initialize turbulent heat fluxes (used with RSM models)
  ! and fluxes coming from interfaces (cases with boiling)
  q_turb(:) = 0.0
  q_int (:) = 0.0

  ! Old values (o and oo)
  if(ini .eq. 1) then
    do c = 1, Grid % n_cells
      t % oo(c) = t % o(c)
      t % o (c) = t % n(c)
    end do
  end if

  ! Gradients
  if(.not. Flow % mass_transfer) then
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

    call Turb_Mod_Face_Cond_And_Stress(turb, con_eff, t_stress, s)
    
    if(Flow % mass_transfer) then
      if(Vof % fun % n(c1) < 0.5 .and.  &
         Vof % fun % n(c2) < 0.5) con_eff = Vof % phase_cond(2)
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

    a12 = a12 - min(v_flux % n(s), 0.0)  &
                 * Flow % capacity(c1) * Flow % density(c1)  ! Flow: 1 -> 2
    a21 = a21 + max(v_flux % n(s), 0.0)  &
                 * Flow % capacity(c2) * Flow % density(c2)  ! Flow: 2 -> 1

    !-----------------------------------------------------!
    !   In case of mass transfer, detach the two phases   !
    !      and add heat transferred to the interface      !
    !-----------------------------------------------------!
    if(Flow % mass_transfer) then
      if(any(Vof % Front % face_at_elem(1:2,s) .ne. 0)) then
        a12  = 0.0
        a21  = 0.0
        f_ex = 0.0  ! included in q_int
        f_im = 0.0  ! phases are detached
        q_int(c1) = q_int(c1) + Vof % q_int(1,s)
        q_int(c2) = q_int(c2) - Vof % q_int(2,s)
      end if
    end if

    ! Cross diffusion part
    t % c(c1) = t % c(c1) + f_ex - f_im
    if(c2 > 0) then
      t % c(c2) = t % c(c2) - f_ex + f_im
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

  !----------------------------------------------!
  !   Explicitly treated diffusion heat fluxes   !
  !   cross diffusion, and heat from interface   !
  !----------------------------------------------!
  do c = 1, Grid % n_cells

    ! Total explicit heat flux
    q_exp = t % c(c) + q_turb(c) + q_int(c)

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
    if(Flow % mass_transfer) then
      if(Vof % fun % n(c) > 0.5) then
        cap_dens(c) = Vof % phase_capa(1) * Vof % phase_dens(1)
      else if(Vof % fun % n(c) < 0.5) then
        cap_dens(c) = Vof % phase_capa(2) * Vof % phase_dens(2)
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

  call Cpu_Timer % Start('Linear_Solver_For_Energy')

  ! Call linear solver to solve the equations
  call Sol % Run(t % solver, t % prec,  &
                 A,                     &
                 t % n,                 &
                 b,                     &
                 t % mniter,            &
                 t % eniter,            &
                 t % tol,               &
                 t % res)

  call Cpu_Timer % Stop('Linear_Solver_For_Energy')

  ! Print some info on the screen
  call Info_Mod_Iter_Fill_At(1, 6, t % name, t % eniter, t % res)

  ! Gradients
  if(.not. Flow % mass_transfer) then
    call Flow % Grad_Variable(t)

  ! If mass transfer, compute gradients with saturation temperature
  else
    call Vof % Grad_Variable_With_Front(t, Vof % t_sat)
    call Flow % Calculate_Grad_Matrix()
  end if

  ! User function
  call User_Mod_End_Of_Compute_Energy(Flow, turb, Vof, Sol, curr_dt, ini)

  call Cpu_Timer % Stop('Compute_Energy (without solvers)')

  end subroutine
