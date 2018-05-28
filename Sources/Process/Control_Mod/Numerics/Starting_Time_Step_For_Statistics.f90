!==============================================================================!
  subroutine Control_Mod_Starting_Time_Step_For_Statistics(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer           :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Int_Item('STARTING_TIME_STEP_FOR_STATISTICS', 1000000, &
                                  val, verbose)

  end subroutine
