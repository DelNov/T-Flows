!==============================================================================!
  subroutine Control_Mod_Number_Of_Phases(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer, intent(out) :: val
  logical, optional    :: verbose
!==============================================================================!

  call Control % Read_Int_Item('NUMBER_OF_PHASES', 1, val, verbose)

  end subroutine
