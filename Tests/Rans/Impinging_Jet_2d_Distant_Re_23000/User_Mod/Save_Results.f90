include '../User_Mod/Impinging_Jet_Nu.f90'
include '../User_Mod/Impinging_Jet_Profiles.f90'

!==============================================================================!
  subroutine User_Mod_Save_Results(grid, save_name) 
!------------------------------------------------------------------------------!
!   Calls User_Impinging_Jet_Nu and User_Impinging_Jet_Profile functions.      !
!------------------------------------------------------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)  :: grid
  character(len=*) :: save_name
!==============================================================================!

  call User_Mod_Impinging_Jet_Nu      (grid, save_name)
  call User_Mod_Impinging_Jet_Profiles(grid, save_name)

  end subroutine
