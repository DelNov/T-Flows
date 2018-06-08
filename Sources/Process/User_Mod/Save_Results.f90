!==============================================================================!
  subroutine User_Mod_Save_Results(grid, name_save)
!------------------------------------------------------------------------------!
!   This is a prototype of customized saving functions.  All cases which       !
!   require special procedures when saving results (1d or 2d profiles, non-    !
!   standard save formats, and alike) should be done here, case-based.         !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
  use Var_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)  :: grid
  character(len=*) :: name_save
!==============================================================================!

  !--------------------------! 
  !   First level comments   !
  !--------------------------! 

  !---------------------------! 
  !   Second level comments   !
  !---------------------------! 

  ! Third level comments

  end subroutine  ! fourth level comments
