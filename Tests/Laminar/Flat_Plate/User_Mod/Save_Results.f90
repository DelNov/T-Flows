include '../User_Mod/Plain_Profiles.f90'
include '../User_Mod/Plain_Nu.f90'

!==============================================================================!
  subroutine User_Mod_Save_Results(Flow, Turb, Vof, Swarm, ts)
!------------------------------------------------------------------------------!
!   Calls user-define subroutines                                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: Flow
  type(Turb_Type),  target :: Turb
  type(Vof_Type),   target :: Vof
  type(Swarm_Type), target :: Swarm
  integer, intent(in)      :: ts
!==============================================================================!

  ! Don't save if this is intial condition, nothing is developed yet
  if(ts .eq. 0) return

  call User_Mod_Plain_Profiles(Flow, Turb, ts)
  call User_Mod_Plain_Nu      (Flow, Turb, ts)

  end subroutine
