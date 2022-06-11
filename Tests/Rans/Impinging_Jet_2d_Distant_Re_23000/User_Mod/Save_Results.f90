include '../User_Mod/Impinging_Jet_Nu.f90'
include '../User_Mod/Impinging_Jet_Profiles.f90'

!==============================================================================!
  subroutine User_Mod_Save_Results(Flow, Turb, Vof, Swarm, ts)
!------------------------------------------------------------------------------!
!   Calls User_Impinging_Jet_Nu and User_Impinging_Jet_Profile functions.      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type)    :: Flow
  type(Turb_Type)     :: Turb
  type(Vof_Type)      :: Vof
  type(Swarm_Type)    :: Swarm
  integer, intent(in) :: ts     ! time step
!==============================================================================!

  ! Don't save if this is intial condition, nothing is developed yet
  if(ts .eq. 0) return

  call User_Mod_Impinging_Jet_Nu      (Turb, ts)
  call User_Mod_Impinging_Jet_Profiles(Turb, ts)

  end subroutine
