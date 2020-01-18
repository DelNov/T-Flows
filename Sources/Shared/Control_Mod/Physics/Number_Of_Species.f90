!==============================================================================!
  subroutine Control_Mod_Number_Of_Species(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer           :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Int_Item('NUMBER_OF_SPECIES', 0, &
                                  val, verbose)

  end subroutine
