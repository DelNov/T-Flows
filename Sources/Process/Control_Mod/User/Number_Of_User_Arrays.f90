!==============================================================================!
  subroutine Control_Mod_Number_Of_User_Arrays(val, verbose)
!------------------------------------------------------------------------------!
!   Reading stuff related to user scalars                                      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer           :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Int_Item('NUMBER_OF_USER_ARRAYS', 0, &
                                  val, verbose)

  end subroutine
