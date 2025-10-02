!==============================================================================!
  subroutine Sgs_Discrete_Random_Walk(Swarm, k, rx, ry, rz)
!------------------------------------------------------------------------------!
!>  Implements the Subgrid Scale (SGS) Discrete Random Walk (DRW) model for
!>  Lagrangian particle tracking in Large Eddy Simulations (LES). This method
!>  introduces stochasticity to the particle motion, reflecting the unresolved
!>  turbulent eddies' effects on particle trajectories.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Eddy lifetime and length: Computes eddy lifetimes and lengths based on   !
!     turbulence parameters, essential for modeling the particle's interaction !
!     with the turbulent flow.                                                 !
!   * Particle time scale: Determines the particle's response time to flow     !
!     changes, crucial for accurately capturing particle dynamics.             !
!   * Random number generation: Utilizes Gaussian random numbers to simulate   !
!     the stochastic nature of turbulence, introducing randomness in particle  !
!     motion.                                                                  !
!   * Turbulent velocity contribution: Calculates turbulent velocity           !
!     fluctuations impacting the particle, combining model constants and       !
!     fluid properties.                                                        !
!                                                                              !
!   Note                                                                       !
!                                                                              !
!   * SGS model introducing stochasticity to LES/k-eps-zeta-f modeled          !
!     quantities seen by particles (from Z. Cheng et. al., (2018) 435-451).    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Swarm_Type), target :: Swarm
  integer,       intent(in) :: k
  real,          intent(in) :: rx
  real,          intent(in) :: ry
  real,          intent(in) :: rz
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),    pointer :: Flow
  type(Grid_Type),     pointer :: Grid
  type(Turb_Type),     pointer :: Turb
  type(Particle_Type), pointer :: Part
  integer,             pointer :: time
  integer :: c                          ! nearest cell
  logical :: flag1, flag2, flag3        ! flags for sigma
  real    :: cd                         ! drag coefficient
  real    :: r1, r2, r3, r4             ! Gaussian randdom numbers
  real    :: l1, l2, l3                 ! length scales
  real    :: zeta1, zeta2, zeta3        ! random numbers (0,1)
  real    :: sigma1, sigma2, sigma3     ! Gaussian random numbers
  real    :: c_o, c_mu                  ! model constants
  real    :: visc_fluid, dens_fluid     ! C/c viscosity
  real    :: t_i_x, t_i_y, t_i_z        ! eddy life time
  real    :: t_e_x, t_e_y, t_e_z        ! instantaneoys eddy life time
  real    :: l_e_x, l_e_y, l_e_z        ! turbulent eddy length
  real    :: t_c_x, t_c_y, t_c_z        ! eddy crossing time for part
  real    :: t_int_x, t_int_y, t_int_z  ! intervals for changing sigma
  real    :: tau_p_x, tau_p_y, tau_p_z  ! particle time scale
!==============================================================================!

  ! Take aliases for Flow
  Flow => Swarm % pnt_flow
  Grid => Swarm % pnt_grid
  Turb => Swarm % pnt_turb
  time => Swarm % time_eim
  Part => Swarm % Particle(k)

  ! Model constants
  c_o   = 3.0  ! calibration constant taken as the same value for k-w model
  c_mu  = 0.09
  r1    = 1.0  ! no seeding (same random numbers will be generated every time!)
  r2    = 1.0
  r3    = 1.0
  flag1 = .false.
  flag2 = .false.
  flag3 = .false.
  t_c_x = 0.0
  t_c_y = 0.0
  t_c_z = 0.0

  ! Nearest cell center to particle
  c = Part % cell

  ! Characteristic velocity and density (needs to be discussed)
  visc_fluid = Flow % viscosity(c)
  dens_fluid = Flow % density(c)

  ! Store it for saving
  Part % dens_fluid = dens_fluid

  ! Particle relaxation time (in this swarm)
  ! this needs to be calculated for the Grid once/ts and not for all
  ! ...particles
  Swarm % tau = Swarm % density * (Swarm % diameter **2) / 18.0 / visc_fluid

  ! Mean life time of the turbulent eddy
  t_i_x = 1.0 / (6.0 * c_mu * abs(Turb % eps % x(c)))
  t_i_y = 1.0 / (6.0 * c_mu * abs(Turb % eps % y(c)))
  t_i_z = 1.0 / (6.0 * c_mu * abs(Turb % eps % z(c)))

  ! Instantaneous eddy lifetime
  call random_number(r1)
  call random_number(r2)
  call random_number(r3)
  call random_number(r4)
  zeta1 = r1  ! random number from 0 to 1
  zeta2 = r2  ! random number from 0 to 1
  zeta3 = r3  ! random number from 0 to 1

  t_e_x = - c_o * log(1.0 - zeta1) * t_i_x
  t_e_y = - c_o * log(1.0 - zeta2) * t_i_y
  t_e_z = - c_o * log(1.0 - zeta3) * t_i_z

  ! Eddy length
  l_e_x = t_e_x * sqrt(TWO_THIRDS * abs(Turb % kin % x(c) * rx))
  l_e_y = t_e_y * sqrt(TWO_THIRDS * abs(Turb % kin % y(c) * ry))
  l_e_z = t_e_z * sqrt(TWO_THIRDS * abs(Turb % kin % z(c) * rz))

  ! Compute drag coefficient
  if (Part % re .ge. 1000.0) then
    cd = 0.43
  else
    cd = 24.0 / Part % re * (Part % f)
  end if

  ! Particle time scale
  tau_p_x = (4.0 * Swarm % density * Swarm % diameter)      &
          / (3.0 * dens_fluid * Part % rel_u_mod * cd)
  tau_p_y = (4.0 * Swarm % density * Swarm % diameter)      &
          / (3.0 * dens_fluid * Part % rel_v_mod * cd)
  tau_p_z = (4.0 * Swarm % density * Swarm % diameter)      &
          / (3.0 * dens_fluid * Part % rel_w_mod * cd)

  l1 = tau_p_x * Part % rel_u_mod
  l2 = tau_p_y * Part % rel_v_mod
  l3 = tau_p_z * Part % rel_w_mod

  ! Eddy crossing time for a particle
  if(l_e_x .lt. l1) then
    t_c_x = - tau_p_x * log(1.0-(l_e_x/Part % rel_u_mod/tau_p_x))
  else if(l_e_y .lt. l2) then
    t_c_y = - tau_p_y * log(1.0-(l_e_y/Part % rel_v_mod/tau_p_y))
  else if(l_e_z .lt. l3) then
    t_c_z = - tau_p_z * log(1.0-(l_e_z/Part % rel_w_mod/tau_p_z))
  end if

  ! Random number interval (for varying sigma)
  if(t_c_x .ne. 0.0) then
    t_int_x = min(t_e_x, t_c_x)
  else
    t_int_x = t_e_x
  end if
  if(t_c_y .ne. 0.0) then
    t_int_y = min(t_e_y, t_c_y)
  else
    t_int_y = t_e_y
  end if
  if(t_c_z .ne. 0.0) then
    t_int_z = min(t_e_z, t_c_z)
  else
    t_int_z = t_e_z
  end if

  ! This can be done in a diff. way ((current_ts - backup_ts) * dt)
  ! Flow time (to check the condition of sigma)
  time = time * Flow % dt

  ! Gaussian random numbers with zero mean and unit std. deviation ...
  ! .. Box-Muller Wiener algorithm was used for this distribution
  ! Re-evaluating sigmas (only) every eddy-cycle interval
  sigma1 = 0.0
  sigma2 = 0.0
  sigma3 = 0.0
  if(r1 .ge. TINY .and. r2 .ge. TINY .and. r3 .ge. TINY) then
    if(time .ge. t_i_x) then
      sigma1 = sqrt(-2.0 * log(r1)) * cos(2.0 * PI * r4)
      flag1  = .true.
    else if(time .ge. t_i_y) then
      sigma2 = sqrt(-2.0 * log(r2)) * cos(2.0 * PI * r4)
      flag2  = .true.
    else if(time .ge. t_i_z) then
      sigma3 = sqrt(-2.0 * log(r3)) * cos(2.0 * PI * r4)
      flag3  = .true.
    end if
  else
    if(First_Proc()) then
      print *, "# Sigma is not continuous!"
    end if
    stop
  end if

  ! Re-initializing the time interval again
  if(flag1 .or. flag2 .or. flag3) then
    Swarm % time_eim = 0
  end if

  ! EIM contribution for the modeled turbulent quantitis
  Part % u_drw = sqrt(TWO_THIRDS * abs(Turb % kin % x(c) * rx)) * sigma1
  Part % v_drw = sqrt(TWO_THIRDS * abs(Turb % kin % y(c) * ry)) * sigma2
  Part % w_drw = sqrt(TWO_THIRDS * abs(Turb % kin % z(c) * rz)) * sigma3

  end subroutine
