!==============================================================================!
  subroutine Swarm_Mod_Sgs_Fukagata(swarm)
!------------------------------------------------------------------------------!
!   SGS model accounting for Brownian diffusion force by Fukagata et al., 2004 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Swarm_Type), target :: swarm
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: flow
  type(Grid_Type),  pointer :: grid
  type(Turb_Type),  pointer :: turb
  type(Var_Type),   pointer :: u, v, w
  integer :: c, c2, l                ! nearest cell
  real    :: r1, r2, zeta            ! randdom number
  real    :: sigma                   ! part vel. std. dev.
  real    :: f_fuka                  ! Fukagata Brownian force
  real    :: lambda, alpha, theta    ! parameters for calc. of sigma
  real    :: c_o, c_eps              ! model constants
  real    :: k_sgs                   ! SGS TKE
  real    :: eps_sgs                 ! SGS dissipation rate
  real    :: t_sgs                   ! SGS fluid lag. integral timescale
  real    :: lf                      ! dynamic model length scale
  real    :: visc_const              ! C/c viscosity
!==============================================================================!

  ! Take aliases for flow
  flow => swarm % pnt_flow
  grid => swarm % pnt_grid
  turb => swarm % pnt_turb
  u    => flow % u
  v    => flow % v
  w    => flow % w

  ! Characteristic viscosity
  ! (You can read it from the control file:
  ! call Control_Mod_Dynamic_Viscosity   (visc_const)
  ! The way it is implemented now it could be different in every processor)
  visc_const = maxval(flow % viscosity(:))

  ! Particle relaxation time (in this swarm)
  ! it only needs to be calculated for the grid once/ts and not for all
  ! ...particles
  swarm % tau = swarm % density * (swarm % diameter **2) / 18.0 / visc_const

  ! LES dynamic SGS TKE and dissipation rate
  do l = 1, grid % n_cells

    if(turb % c_dyn(l) .eq. 0.0) then
      swarm % f_fuka_x(l) = 0.0
      swarm % f_fuka_y(l) = 0.0
      swarm % f_fuka_z(l) = 0.0

    else ! Grid is allowing the model to add a contribution 

      ! Length scale of Dynamic model 
      lf = grid % vol(l) ** ONE_THIRD

      ! SGS turbulent kinetic energy
      k_sgs = lf*lf                &  ! delta^2
            * turb % c_dyn(l)      &  ! c_dynamic
            * (flow % shear(l))    &  ! |S|^2
            * (flow % shear(l))


      ! SGS dissipation rate
      eps_sgs = c_eps * (k_sgs**1.5) / lf

      ! SGS fluid lagrangian integral timescale ...
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
        print *, "zeta is NOT continuous!"
        stop
      end if

      ! Brownian diffusion force following isotropic turbulence
      swarm % f_fuka_x(l) = sigma * zeta / swarm % dt
      swarm % f_fuka_y(l) = sigma * zeta / swarm % dt
      swarm % f_fuka_z(l) = sigma * zeta / swarm % dt

    end if
  end do

  end subroutine
