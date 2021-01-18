!==============================================================================!
  subroutine Compute_Energy(flow, turb, mult, sol, ini)
!------------------------------------------------------------------------------!
!   Purpose: Solve transport equation for scalar (such as temperature)         !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use User_Mod
  use Work_Mod, only: capacity_x_density => r_cell_11
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
  real                       :: a12, a21, con_eff_f, con_t_f
  real                       :: f_ex, f_im, tx_f, ty_f, tz_f
  real                       :: pr_t1, pr_t2, pr_tf
  real                       :: t_stress, dt
  real                       :: cap_dens_c1, cap_dens_c2
  real                       :: ut_x_cap_dens_s, &
                                vt_x_cap_dens_s, &
                                wt_x_cap_dens_s
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
  call Field_Mod_Alias_Momentum  (flow, u, v, w)
  call Field_Mod_Alias_Energy    (flow, t)
  call Turb_Mod_Alias_Heat_Fluxes(turb, ut, vt, wt)
  call Solver_Mod_Alias_System   (sol, a, b)

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

  !---------------!
  !               !
  !   Advection   !
  !               !
  !---------------!
  call Numerics_Mod_Advection_Term(t, flow % capacity * flow % density,  &
                                   v_flux % n, sol)

  !--------------!
  !              !
  !   Difusion   !
  !              !
  !--------------!

  !----------------------------!
  !   Spatial discretization   !
  !----------------------------!
  if(turb % model .ne. NO_TURBULENCE_MODEL .and.  &
     turb % model .ne. DNS) then
  end if

  do s = 1, grid % n_faces

    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    if(turb % model .ne. LES_SMAGORINSKY     .and.  &
       turb % model .ne. LES_DYNAMIC         .and.  &
       turb % model .ne. HYBRID_LES_PRANDTL  .and.  &
       turb % model .ne. LES_WALE            .and.  &
       turb % model .ne. NO_TURBULENCE_MODEL .and.  &
       turb % model .ne. DNS) then
      pr_t1 = Turb_Mod_Prandtl_Number(turb, c1)
      pr_t2 = Turb_Mod_Prandtl_Number(turb, c2)
      pr_tf = grid % fw(s) * pr_t1 + (1.0-grid % fw(s)) * pr_t2
    else
      pr_tf = pr_t
    end if

    ! Gradients on the cell face (fw corrects situation close to the wall)
    tx_f = grid % fw(s) * t % x(c1) + (1.0-grid % fw(s)) * t % x(c2)
    ty_f = grid % fw(s) * t % y(c1) + (1.0-grid % fw(s)) * t % y(c2)
    tz_f = grid % fw(s) * t % z(c1) + (1.0-grid % fw(s)) * t % z(c2)
    if(turb % model .ne. NO_TURBULENCE_MODEL .and.  &
       turb % model .ne. DNS) then
      con_eff_f =                                                              &
               grid % fw(s) * (flow % conductivity(c1) +                       &
                               flow % capacity(c1) * turb % vis_t(c1) / pr_tf) &
        + (1.0-grid % fw(s))* (flow % conductivity(c2) +                       &
                               flow % capacity(c2) * turb % vis_t(c2) / pr_tf)
      con_t_f  = grid % fw(s) *flow % capacity(c1) * turb % vis_t(c1) / pr_tf  &
          + (1.0-grid % fw(s))*flow % capacity(c2) * turb % vis_t(c2) / pr_tf
    else
      con_eff_f =                                      &
               grid % fw(s) * flow % conductivity(c1)  &
        + (1.0-grid % fw(s))* flow % conductivity(c2)
    end if
    if(turb % model .eq. HYBRID_LES_RANS) then
      con_eff_f =                                                              &
               grid % fw(s) * (flow % conductivity(c1) +                       &
                               flow % capacity(c1) * turb % vis_t_eff(c1)      &
                               / pr_tf)                                        &
       + (1.0-grid % fw(s)) * (flow % conductivity(c2) +                       &
                               flow % capacity(c2) * turb % vis_t_eff(c2)      &
                               / pr_tf)
     con_t_f  = grid % fw(s) * flow % capacity(c1) * turb % vis_t_eff(c1)      &
                               / pr_tf                                         &
        + (1.0-grid % fw(s)) * flow % capacity(c2) * turb % vis_t_eff(c2)      &
                               / pr_tf
    end if

    if(turb % model .eq. K_EPS        .or.  &
       turb % model .eq. K_EPS_ZETA_F .or.  &
       turb % model .eq. HYBRID_LES_RANS) then
      if(c2 < 0) then
        if(Var_Mod_Bnd_Cond_Type(t, c2) .eq. WALL .or.  &
           Var_Mod_Bnd_Cond_Type(t, c2) .eq. WALLFL) then
          con_eff_f = turb % con_w(c1)
        end if
      end if
    end if

    ! Total (exact) diffusion flux
    ! This is last term in equation 2.33 in Denner's thesis because:
    ! grid % sx (T-Flows) == n_f * A_f (Denner)
    f_ex = con_eff_f * (   tx_f * grid % sx(s)   &
                         + ty_f * grid % sy(s)   &
                         + tz_f * grid % sz(s))

    ! Implicit part of the diffusion flux, treated by linear system
    ! This is also term in equation 2.33 in Denner's thesis because:
    ! a % fc * grid % dx (T-Flows) == alpha_f * s_f * A_f (Denner)
    f_im = con_eff_f * a % fc(s)         &
          * (   tx_f * grid % dx(s)      &
              + ty_f * grid % dy(s)      &
              + tz_f * grid % dz(s) )

    ! Cross diffusion part
    t % c(c1) = t % c(c1) + f_ex - f_im
    if(c2 > 0) then
      t % c(c2) = t % c(c2) - f_ex + f_im
    end if

    !---------------------------!
    !                           !
    !   Turbulent heat fluxes   !
    !                           !
    !---------------------------!
    if(turb % model .eq. RSM_HANJALIC_JAKIRLIC .or.  &
       turb % model .eq. RSM_MANCEAU_HANJALIC) then

      cap_dens_c1 = flow % capacity(c1) * flow % density(c1)
      cap_dens_c2 = flow % capacity(c2) * flow % density(c2)

      ! Turbulent heat fluxes according to GGDH scheme
      ! (first line is GGDH, second line is SGDH substratced
      ut_x_cap_dens_s =  (    grid % fw(s)  * ut % n(c1) * cap_dens_c1   &
                      +  (1.0-grid % fw(s)) * ut % n(c2) * cap_dens_c2)
      vt_x_cap_dens_s =  (    grid % fw(s)  * vt % n(c1) * cap_dens_c1   &
                      +  (1.0-grid % fw(s)) * vt % n(c2) * cap_dens_c2)
      wt_x_cap_dens_s =  (    grid % fw(s)  * wt % n(c1) * cap_dens_c1   &
                      +  (1.0-grid % fw(s)) * wt % n(c2) * cap_dens_c2)
      t_stress = - (  ut_x_cap_dens_s * grid % sx(s)                     &
                    + vt_x_cap_dens_s * grid % sy(s)                     &
                    + wt_x_cap_dens_s * grid % sz(s) )                   &
                    - (con_t_f * (  tx_f * grid % sx(s)     &
                                  + ty_f * grid % sy(s)     &
                                  + tz_f * grid % sz(s)) )

      ! Put the influence of turbulent heat fluxes explicitly in the system
      b(c1) = b(c1) + t_stress
      if(c2 > 0) then
        b(c2) = b(c2) - t_stress
      end if
    end if  ! if models are of RSM type

    ! Calculate the coefficients for the sysytem matrix
    a12 = con_eff_f * a % fc(s)
    a21 = con_eff_f * a % fc(s)

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
    capacity_x_density(c) = flow % capacity(c) * flow % density(c)
  end do
  call Numerics_Mod_Inertial_Term(t, capacity_x_density, sol, dt)

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

  if (mult % phase_change) then
    do c = 1, grid % n_cells
      if (mult % ic(c) == 2) then
        a % val(a % dia(c)) = HUGE
        b(c) = HUGE * flow % sat_temperature
        write(*,*) c, a % val(a % dia(c)), b(c), flow % sat_temperature, t % n(c), t % urf
      else
        b(c) = b(c) + mult % qci(c)
      end if
    end do
  end if

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
