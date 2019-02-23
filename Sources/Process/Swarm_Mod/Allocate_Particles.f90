!==============================================================================!
  subroutine Swarm_Mod_Allocate(flow, swarm)
!------------------------------------------------------------------------------!
!   Allocates memory to store the charge of each particle                      !
!   It assumes that the number of particles was read from the control file     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: flow
  type(Swarm_Type), target :: swarm
!==============================================================================!

  ! Take aliases to object particle flow around
  swarm % pnt_flow => flow
  swarm % pnt_grid => flow % pnt_grid

  ! Allocate memory for all of them
  allocate(swarm % particles(swarm % n_particles))

  end subroutine
