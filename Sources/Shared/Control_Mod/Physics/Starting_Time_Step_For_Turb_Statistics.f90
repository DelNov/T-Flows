!==============================================================================!
  subroutine Control_Mod_Starting_Time_Step_For_Turb_Statistics(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer, intent(out) :: val
  logical, optional    :: verbose
!==============================================================================!

  call Control % Read_Int_Item('STARTING_TIME_STEP_FOR_TURB_STATISTICS',  &
                                HUGE_INT, val, verbose)

  end subroutine
