!==============================================================================!
  subroutine Swarm_Mod_Allocate(flow, swarm, n_part)
!------------------------------------------------------------------------------!
!   Allocates memory to store the charge of each particle                      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: flow
  type(Swarm_Type), target :: swarm
  integer                  :: n_part
!-----------------------------------[Locals]-----------------------------------!
!==============================================================================!

  ! Take aliases to object particle flow around
  swarm % pnt_flow => flow
  swarm % pnt_grid => flow % pnt_grid

  ! Store number of particles (to be assigned by user later!)
  swarm % n_particles = n_part

  ! Allocate memory for all of them
  allocate(swarm % particles(n_part))

  end subroutine
