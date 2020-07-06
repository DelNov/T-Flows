!==============================================================================!
  subroutine Swarm_Mod_Sgs_Fukagata(swarm)
!------------------------------------------------------------------------------!
!   SGS model accounting for Brownian diffusion force by Fukagata et al., 2004 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Swarm_Type), target :: swarm
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),    pointer :: flow
  type(Grid_Type),     pointer :: grid
  type(Turb_Type),     pointer :: turb
  type(Var_Type),      pointer :: u, v, w
  integer                      :: c, c2                 ! nearest cell
  real                         :: r1, r2, zeta          ! random number
  real                         :: sigma                 ! part vel. std. dev.
  real                         :: f_fuka                ! Fukagata Brownian f.
  real                         :: lambda, alpha, theta  ! for calc. of sigma
  real                         :: c_o, c_eps            ! model constants
  real                         :: k_sgs                 ! SGS kin
  real                         :: eps_sgs               ! SGS eps
  real                         :: t_sgs                 ! SGS time scale
  real                         :: lf                    ! dynamic model l scale
  real                         :: visc_const            ! const viscosity
!==============================================================================!

  ! Take aliases for flow
  flow => swarm % pnt_flow
  grid => swarm % pnt_grid
  turb => swarm % pnt_turb
  u    => flow % u
  v    => flow % v
  w    => flow % w

  ! Model constants
  c_o   = 2.1
  c_eps = 1.0
  r1    = 1.0
  r2    = 1.0

  ! Update eddy viscosity from the dynamic model
  call Turb_Mod_Vis_T_Dynamic(turb)

  ! Characteristic viscosity
  visc_const = maxval(flow % viscosity(:))

  ! Particle relaxation time (in this swarm)
  ! it only needs to be calculated for the grid once/ts and not for all ...
  ! ...particles
  swarm % tau = swarm % density * (swarm % diameter **2) / 18.0 / visc_const

  ! LES dynamic SGS TKE and dissipation rate
  do c = 1, grid % n_cells

    if(turb % c_dyn(c) .eq. 0.0) then
      swarm % f_fuka_x(c) = 0.0
      swarm % f_fuka_y(c) = 0.0
      swarm % f_fuka_z(c) = 0.0

    else ! Grid is allowing the model to add a contribution

      ! Length scale of Dynamic model
      lf = grid % vol(c) ** ONE_THIRD

      ! SGS turbulent kinetic energy
      k_sgs = lf*lf                &  ! delta^2
            * turb % c_dyn(c)      &  ! c_dynamic
            * (flow % shear(c))    &  ! |S|^2
            * (flow % shear(c))


      ! SGS dissipation rate
      eps_sgs = c_eps * (k_sgs**1.5) / lf

      ! SGS fluid Lagrangian integral timescale ...
      ! ... along an inertial particle's path
      t_sgs = (1.0/(0.5 + 0.75*c_o))*(k_sgs/eps_sgs)

      ! Standard deviation of particle velocity due to SGS vel. fluctuations 
      theta  = swarm % tau / t_sgs
      alpha  = flow % dt / swarm % tau
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
      swarm % f_fuka_x(c) = sigma * zeta / swarm % dt
      swarm % f_fuka_y(c) = sigma * zeta / swarm % dt
      swarm % f_fuka_z(c) = sigma * zeta / swarm % dt

    end if
  end do

  end subroutine
