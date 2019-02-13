include '../User_Mod/Impinging_Jet_Nu.f90'
include '../User_Mod/Impinging_Jet_Profiles.f90'

!==============================================================================!
  subroutine User_Mod_Save_Results(flow, save_name)
!------------------------------------------------------------------------------!
!   Calls User_Impinging_Jet_Nu and User_Impinging_Jet_Profile functions.      !
!------------------------------------------------------------------------------!
  use Field_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type) :: flow
  character(len=*) :: save_name
!==============================================================================!

  call User_Mod_Impinging_Jet_Nu      (flow, save_name)
  call User_Mod_Impinging_Jet_Profiles(flow, save_name)

  end subroutine
