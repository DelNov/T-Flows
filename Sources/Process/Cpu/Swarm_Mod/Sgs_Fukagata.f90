!==============================================================================!
  subroutine Sgs_Fukagata(Swarm)
!------------------------------------------------------------------------------!
!> Implements the Subgrid Scale (SGS) Brownian diffusion force model by
!> Fukagata et al., 2004, for Lagrangian particle tracking in turbulence
!> simulations. This model is particularly effective for representing the
!> impacts of unresolved turbulent eddies on particle behavior in LES.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Turbulence modeling: Applies the SGS model to calculate turbulence       !
!     parameters like kinetic energy and dissipation rate, essential for       !
!     accurate turbulence representation.                                      !
!   * Particle relaxation time: Computes the particle's response time scale    !
!     to flow changes, crucial for capturing particle dynamics within the      !
!     turbulent flow.                                                          !
!   * Brownian diffusion force: Calculates Brownian diffusion force components !
!     impacting particle motion, accounting for turbulence-induced randomness. !
!   * Random number generation: Utilizes Gaussian random numbers to model the  !
!     stochastic nature of turbulence, impacting particle trajectories.        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Swarm_Type), target :: Swarm  !! the swarm of particles
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: Flow
  type(Grid_Type),  pointer :: Grid
  type(Turb_Type),  pointer :: Turb
  type(Var_Type),   pointer :: u, v, w
  integer                   :: c                     ! nearest cell
  real                      :: r1, r2, zeta          ! random number
  real                      :: sigma                 ! part vel. std. dev.
  real                      :: lambda, alpha, theta  ! for calc. of sigma
  real                      :: c_o, c_eps            ! model constants
  real                      :: k_sgs                 ! SGS kin
  real                      :: eps_sgs               ! SGS eps
  real                      :: t_sgs                 ! SGS time scale
  real                      :: lf                    ! dynamic model l scale
  real                      :: visc_const            ! const viscosity
!==============================================================================!

  ! Take aliases for Flow
  Flow => Swarm % pnt_flow
  Grid => Swarm % pnt_grid
  Turb => Swarm % pnt_turb
  u    => Flow % u
  v    => Flow % v
  w    => Flow % w

  ! Model constants
  c_o   = 2.1
  c_eps = 1.0
  r1    = 1.0
  r2    = 1.0

  ! Update eddy viscosity from the dynamic model
  call Turb % Vis_T_Dynamic()

  ! Characteristic viscosity
  visc_const = maxval(Flow % viscosity(:))

  ! Particle relaxation time (in this Swarm)
  ! it only needs to be calculated for the Grid once/ts and not for all ...
  ! ...particles
  Swarm % tau = Swarm % density * (Swarm % diameter **2) / 18.0 / visc_const

  ! LES dynamic SGS TKE and dissipation rate
  do c = Cells_In_Domain_And_Buffers()

    if(Turb % c_dyn(c) .eq. 0.0) then
      Swarm % f_fuka_x(c) = 0.0
      Swarm % f_fuka_y(c) = 0.0
      Swarm % f_fuka_z(c) = 0.0

    else ! Grid is allowing the model to add a contribution

      ! Length scale of Dynamic model
      lf = Grid % vol(c) ** ONE_THIRD

      ! SGS turbulent kinetic energy
      k_sgs = lf*lf                &  ! delta^2
            * Turb % c_dyn(c)      &  ! c_dynamic
            * (Flow % shear(c))    &  ! |S|^2
            * (Flow % shear(c))


      ! SGS dissipation rate
      eps_sgs = c_eps * (k_sgs**1.5) / lf

      ! SGS fluid Lagrangian integral timescale ...
      ! ... along an inertial particle's path
      t_sgs = (1.0/(0.5 + 0.75*c_o))*(k_sgs/eps_sgs)

      ! Standard deviation of particle velocity due to SGS vel. fluctuations 
      theta  = Swarm % tau / t_sgs
      alpha  = Flow % dt / Swarm % tau
      lambda = (1.0/1.0 + theta)                    &
             * (1.0-exp(-1.0*alpha*(1.0+theta)))    &
             - (1.0/1.0 - theta) * exp(-2.0*alpha)  &
             * (1.0-exp(alpha*(1.0-theta)))

      ! Velocity standard deviation
      sigma  = sqrt(TWO_THIRDS * k_sgs * lambda)

      ! Gaussian random number with zero mean and unit variance
      call random_number(r1)
      call random_number(r2)

      ! Box-Muller Wiener (Gaussian) random distribution
      if(r1 .ge. TINY) then
        zeta = sqrt(-2.0 * log(r1)) * cos(2.0 * PI * r2)
      else
        print *, "# zeta is not continuous!"
        print *, "# This error is critical; exiting!"
        stop
      end if

      ! Brownian diffusion force following isotropic turbulence
      Swarm % f_fuka_x(c) = sigma * zeta / Swarm % dt
      Swarm % f_fuka_y(c) = sigma * zeta / Swarm % dt
      Swarm % f_fuka_z(c) = sigma * zeta / Swarm % dt

    end if
  end do

  end subroutine
