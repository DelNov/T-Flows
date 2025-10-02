!==============================================================================!
  subroutine User_Mod_Beginning_Of_Simulation(Flow, Turb, Vof, Swarm)
!------------------------------------------------------------------------------!
!   This function is called at the beginning of simulation.                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: Flow
  type(Turb_Type),  target :: Turb
  type(Vof_Type),   target :: Vof
  type(Swarm_Type), target :: Swarm
!==============================================================================!

  end subroutine
