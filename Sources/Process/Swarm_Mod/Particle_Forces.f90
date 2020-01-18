!==============================================================================!
  subroutine Swarm_Mod_Particle_Forces(swarm, k)
!------------------------------------------------------------------------------!
!                 Computes the forces exerted on the particle                  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Swarm_Type), target :: swarm
  integer                  :: k     ! particle number
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),    pointer :: flow
  type(Particle_Type), pointer :: part
  real                         :: cd          ! drag coefficient
  real                         :: part_vol    ! particle volume
  real                         :: part_area   ! particle area
  real                         :: dens_const  ! characteristic density 
!==============================================================================!

  ! Take aliases
  flow => swarm % pnt_flow
  part => swarm % particle(k)

  ! Particle surface area (assuming spherical shape)
  part_area =  PI * (part % d ** 2)

  ! Particle volume
  part_vol = ONE_SIXTH * PI * (part % d ** 3)

  !----------------------------------------------------!
  !  Compute the buoyancy force (acting in +ve y-dir)  !
  !----------------------------------------------------!

  ! if the particle is lighter than the fluid (i.e. fb is +ve)...
  ! ...then both buoyancy and drag forces will act on the particle...
  ! ...in the same direction (upwards), and in this case the particle...
  ! ...will be moving upwards as well. Otherwise fb will be -ve and it...
  ! ...will be deducted from the total force.

  ! Characteristic density (needs to be discussed):
  dens_const = maxval(flow % density(:))

  part % fb_x = 0.0
  part % fb_y = (dens_const - part % density) * EARTH_G * part_vol
  part % fb_y = 0.0

  ! Compute drag coefficient
  if (part % re .ge. 1000.0) then
    cd = 0.43
  else
    cd = 24.0 / part % re * (part % f)
  end if

  !-------------------------------------------------------------------!
  !   Compute the drag force (acting in particle counter direction)   !
  !-------------------------------------------------------------------!
  part % fd_x = .5 * cd * dens_const * part_area * part % rel_vel * part % rel_u
  part % fd_y = .5 * cd * dens_const * part_area * part % rel_vel * part % rel_v
  part % fd_z = .5 * cd * dens_const * part_area * part % rel_vel * part % rel_w

  !-----------------------------!
  !   Compute the total force   !
  !-----------------------------!
  part % ft_x = part % fb_x + part % fd_x
  part % ft_y = part % fb_y + part % fd_y
  part % ft_z = part % fb_z + part % fd_z

  end subroutine
