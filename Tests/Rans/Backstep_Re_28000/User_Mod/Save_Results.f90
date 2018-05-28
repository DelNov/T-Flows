include '../User_Mod/Backstep_Profiles.f90'
include '../User_Mod/Backstep_Cf_St.f90'

!==============================================================================!
  subroutine User_Mod_Save_Results(grid, save_name) 
!------------------------------------------------------------------------------!
!   Calls User_Backstep_Profiles and User_Backstep_Cf_St functions.            !
!------------------------------------------------------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)  :: grid
  character(len=*) :: save_name
!==============================================================================!

  call User_Mod_Backstep_Profiles(grid, save_name)
  call User_Mod_Backstep_Cf_St   (grid, save_name)

  end subroutine
