!==============================================================================!
  subroutine User_Mod_Beginning_Of_Time_Step(flow, turb, swarm, n, time)
!------------------------------------------------------------------------------!
!   This function is called at the beginning of time step.                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: flow
  type(Turb_Type),  target :: turb
  type(Swarm_Type), target :: swarm
  integer                  :: n     ! time step
  real                     :: time  ! physical time
!==============================================================================!

  end subroutine
