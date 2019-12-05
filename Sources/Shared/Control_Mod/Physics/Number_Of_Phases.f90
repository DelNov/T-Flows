!==============================================================================!
  subroutine Control_Mod_Number_Of_Phases(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer           :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Int_Item('NUMBER_OF_PHASES', 1, &
                                  val, verbose)

  end subroutine
