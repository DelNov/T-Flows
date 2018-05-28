!==============================================================================!
  subroutine Control_Mod_Results_Save_Interval(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer           :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Int_Item('RESULTS_SAVE_INTERVAL', 50, &
                                  val, verbose)

  end subroutine
