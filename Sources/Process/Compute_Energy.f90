!==============================================================================!
  subroutine Compute_Energy(flow, turb, mult, sol, ini)
!------------------------------------------------------------------------------!
!   Purpose: Solve transport equation for scalar (such as temperature)         !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use User_Mod
  use Work_Mod, only: cap_dens => r_cell_11
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Arguments]--------------------------------!
  type(Field_Type),      target :: flow
  type(Turb_Type),       target :: turb
  type(Multiphase_Type), target :: mult
  type(Solver_Type),     target :: sol
  integer                       :: ini
!-----------------------------------[Locals]-----------------------------------! 
  type(Grid_Type),   pointer :: grid
  type(Var_Type),    pointer :: u, v, w, t
  type(Var_Type),    pointer :: ut, vt, wt
  type(Face_Type),   pointer :: v_flux
  type(Matrix_Type), pointer :: a
  real, contiguous,  pointer :: b(:)
  integer                    :: c, s, c1, c2
  real                       :: a12, a21, con_eff
  real                       :: f_ex, f_im, tx_f, ty_f, tz_f
  real                       :: t_stress, dt
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
!   left  hand s.               a                 [J/(s K)]
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
!   heat flux                   flow % heat_flux  [W/m^2]
!==============================================================================!

  call Cpu_Timer_Mod_Start('Compute_Energy (without solvers)')

  ! Take aliases
  grid   => flow % pnt_grid
  v_flux => flow % v_flux
  dt     =  flow % dt
  call Field_Mod_Alias_Momentum(flow, u, v, w)
  call Field_Mod_Alias_Energy  (flow, t)
  call Solver_Mod_Alias_System (sol, a, b)

  ! User function
  call User_Mod_Beginning_Of_Compute_Energy(flow, turb, mult, ini)

  ! Initialize matrix and right hand side
  a % val(:) = 0.0
  b      (:) = 0.0

  ! Old values (o and oo)
  if(ini .eq. 1) then
    do c = 1, grid % n_cells
      t % oo(c) = t % o(c)
      t % o (c) = t % n(c)
    end do
  end if

  ! Gradients
  call Field_Mod_Grad_Variable(flow, t)

  ! Compute helping variable
  do c = -grid % n_bnd_cells, grid % n_cells
    cap_dens(c) = flow % capacity(c) * flow % density(c)
  end do

  !---------------!
  !               !
  !   Advection   !
  !               !
  !---------------!
  call Numerics_Mod_Advection_Term(t, cap_dens, v_flux % n, sol)

  !--------------!
  !              !
  !   Difusion   !
  !              !
  !--------------!

  !----------------------------!
  !   Spatial discretization   !
  !----------------------------!
  do s = 1, grid % n_faces

    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    call Turb_Mod_Face_Cond_And_Stress(turb, con_eff, t_stress, s)

    ! Gradients on the cell face (fw corrects situation close to the wall)
    tx_f = grid % fw(s) * t % x(c1) + (1.0-grid % fw(s)) * t % x(c2)
    ty_f = grid % fw(s) * t % y(c1) + (1.0-grid % fw(s)) * t % y(c2)
    tz_f = grid % fw(s) * t % z(c1) + (1.0-grid % fw(s)) * t % z(c2)

    ! Total (exact) diffusion flux
    ! This is last term in equation 2.33 in Denner's thesis because:
    ! grid % sx (T-Flows) == n_f * A_f (Denner)
    f_ex = con_eff * (   tx_f * grid % sx(s)   &
                       + ty_f * grid % sy(s)   &
                       + tz_f * grid % sz(s))

    ! Implicit part of the diffusion flux, treated by linear system
    ! This is also term in equation 2.33 in Denner's thesis because:
    ! a % fc * grid % dx (T-Flows) == alpha_f * s_f * A_f (Denner)
    f_im = con_eff * a % fc(s) * (   tx_f * grid % dx(s)    &
                                   + ty_f * grid % dy(s)    &
                                   + tz_f * grid % dz(s) )

    ! Cross diffusion part
    t % c(c1) = t % c(c1) + f_ex - f_im
    if(c2 > 0) then
      t % c(c2) = t % c(c2) - f_ex + f_im
    end if

    ! Put the influence of turbulent heat fluxes explicitly in the system
    b(c1) = b(c1) + t_stress
    if(c2 > 0) then
      b(c2) = b(c2) - t_stress
    end if

    !------------------------------------------------------!
    !                                                      !
    !   Calculate the coefficients for the system matrix   !
    !                                                      !
    !------------------------------------------------------!
    a12 = con_eff * a % fc(s)
    a21 = con_eff * a % fc(s)

    a12 = a12 - min(v_flux % n(s), 0.0)  &
                 * flow % capacity(c1) * flow % density(c1)  ! flow: 1 -> 2
    a21 = a21 + max(v_flux % n(s), 0.0)  &
                 * flow % capacity(c2) * flow % density(c2)  ! flow: 2 -> 1

    ! Fill the system matrix
    if(c2 > 0) then
      a % val(a % dia(c1))  = a % val(a % dia(c1)) + a12
      a % val(a % dia(c2))  = a % val(a % dia(c2)) + a21
      a % val(a % pos(1,s)) = a % val(a % pos(1,s)) - a12
      a % val(a % pos(2,s)) = a % val(a % pos(2,s)) - a21
    else if(c2 .lt. 0) then
      ! Outflow is included because of the flux
      ! corrections which also affects velocities
      if( (Var_Mod_Bnd_Cond_Type(t, c2) .eq. INFLOW) .or.  &
          (Var_Mod_Bnd_Cond_Type(t, c2) .eq. WALL)   .or.  &
          (Var_Mod_Bnd_Cond_Type(t, c2) .eq. CONVECT) ) then
        a % val(a % dia(c1)) = a % val(a % dia(c1)) + a12
        b(c1)  = b(c1)  + a12 * t % n(c2)
      ! In case of wallflux 
      else if(Var_Mod_Bnd_Cond_Type(t, c2) .eq. WALLFL) then
        b(c1) = b(c1) + grid % s(s) * t % q(c2)
      end if
    end if

  end do  ! through sides

  ! Cross diffusion terms are treated explicity
  do c = 1, grid % n_cells
    if(t % c(c) >= 0) then
      b(c)  = b(c) + t % c(c)
    else
      a % val(a % dia(c)) = a % val(a % dia(c))  &
                          - t % c(c) / (t % n(c) + MICRO)
    end if
  end do

  !--------------------!
  !                    !
  !   Inertial terms   !
  !                    !
  !--------------------!
  do c = -grid % n_bnd_cells, grid % n_cells
    cap_dens(c) = flow % capacity(c) * flow % density(c)
  end do
  call Numerics_Mod_Inertial_Term(t, cap_dens, sol, dt)

  !--------------------!
  !                    !
  !   User source(s)   !
  !                    !
  !--------------------!
  call User_Mod_Source(flow, t, a, b)

  !-------------------------------!
  !                               !
  !   Solve the equations for t   !
  !                               !
  !-------------------------------!

  ! Under-relax the equations
  call Numerics_Mod_Under_Relax(t, sol)

! if (mult % mass_transfer) then
!   do c = 1, grid % n_cells
!     if (mult % ic(c) == 2) then
!       a % val(a % dia(c)) = HUGE
!       b(c) = HUGE * flow % sat_temperature
!       write(*,*) c, a % val(a % dia(c)), b(c), flow % sat_temperature, t % n(c), t % urf
!     else
!       b(c) = b(c) + mult % qci(c)
!     end if
!   end do
! end if

  ! Call linear solver to solve the equations
  call Cpu_Timer_Mod_Start('Linear_Solver_For_Energy')
  call Solver_Mod_Bicg(sol,          &
                       t % n,        &
                       b,            &
                       t % precond,  &
                       t % mniter,   &
                       t % eniter,   &
                       t % tol,      &
                       t % res)
  call Cpu_Timer_Mod_Stop('Linear_Solver_For_Energy')

  ! Print some info on the screen
  call Info_Mod_Iter_Fill_At(1, 6, t % name, t % eniter, t % res)

  call Field_Mod_Grad_Variable(flow, t)

  ! User function
  call User_Mod_End_Of_Compute_Energy(flow, turb, mult, ini)

  call Cpu_Timer_Mod_Stop('Compute_Energy (without solvers)')

  end subroutine
