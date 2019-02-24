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
  real                         :: cd         ! drag coefficient
  real                         :: part_vol   ! particle volume
  real                         :: part_area  ! particle area
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

  part % fb = (density - part % density) * EARTH_G * part_vol

  ! Compute drag coefficient
  if (part % re .ge. 1000.0) then
    cd = 0.43
  else
    cd = 24.0 / part % re * (part % f)
  end if

  !-------------------------------------------------------------------!
  !   Compute the drag force (acting in particle counter direction)   !
  !-------------------------------------------------------------------!
  part % fd = .5 * cd * density * part_area  &
            * abs(part % rel_vv) * part % rel_vv

  !-----------------------------!
  !   Compute the total force   !
  !-----------------------------!
  part % ft = part % fb + part % fd

  end subroutine
