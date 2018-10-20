!==============================================================================!
  subroutine Control_Mod_Starting_Time_Step_For_Statistics(val, verbose)
!------------------------------------------------------------------------------!
!---------------------------------[Modules]------------------------------------!
  use Const_Mod, only: HUGE_INT
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer           :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Int_Item('STARTING_TIME_STEP_FOR_STATISTICS',  &
                                 HUGE_INT, val, verbose)

  end subroutine
