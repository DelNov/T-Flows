include '../User_Mod/Impinging_Jet_Nu.f90'
include '../User_Mod/Impinging_Jet_Profiles.f90'

!==============================================================================!
  subroutine User_Mod_Save_Results(Flow, turb, Vof, swarm, ts)
!------------------------------------------------------------------------------!
!   Calls User_Impinging_Jet_Nu and User_Impinging_Jet_Profile functions.      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type)    :: Flow
  type(Turb_Type)     :: turb
  type(Vof_Type)      :: Vof
  type(Swarm_Type)    :: swarm
  integer, intent(in) :: ts     ! time step
!==============================================================================!

  call User_Mod_Impinging_Jet_Nu      (turb)
  call User_Mod_Impinging_Jet_Profiles(turb)

  end subroutine
