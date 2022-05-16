!==============================================================================!
  subroutine Particle_Forces(Swarm, k)
!------------------------------------------------------------------------------!
!   Computes the forces exerted on the particle                                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Swarm_Type), target :: Swarm
  integer, intent(in)       :: k     ! particle number
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),    pointer :: Flow
  type(Particle_Type), pointer :: Part
  real                         :: cd          ! drag coefficient
  real                         :: part_vol    ! particle volume
  real                         :: part_area   ! particle area
  real                         :: dens_fluid  ! characteristic density 
!==============================================================================!

  ! Take aliases
  Flow => Swarm % pnt_flow
  Part => Swarm % Particle(k)

  ! Particle surface area (assuming spherical shape)
  part_area =  PI * (Part % d ** 2)

  ! Particle volume
  part_vol = ONE_SIXTH * PI * (Part % d ** 3)

  !----------------------------------------------------!
  !  Compute the buoyancy force (acting in +ve y-dir)  !
  !----------------------------------------------------!

  ! if the particle is lighter than the fluid (i.e. fb is +ve)...
  ! ...then both buoyancy and drag forces will act on the particle...
  ! ...in the same direction (upwards), and in this case the particle...
  ! ...will be moving upwards as well. Otherwise fb will be -ve and it...
  ! ...will be deducted from the total force.

  ! Characteristic density (needs to be discussed):
  dens_fluid = Flow % density(Part % cell)

  ! Store it for future saving
  Part % dens_fluid = dens_fluid

  ! Compute drag coefficient
  if (Part % re .ge. 1000.0) then
    cd = 0.43
  else
    cd = 24.0 / Part % re * (Part % f)
  end if

  !-------------------------------------------------------------------!
  !   Compute the drag force (acting in particle counter direction)   !
  !-------------------------------------------------------------------!
  Part % fd_x = .5 * cd * dens_fluid * part_area * Part % rel_vel * Part % rel_u
  Part % fd_y = .5 * cd * dens_fluid * part_area * Part % rel_vel * Part % rel_v
  Part % fd_z = .5 * cd * dens_fluid * part_area * Part % rel_vel * Part % rel_w

  end subroutine
