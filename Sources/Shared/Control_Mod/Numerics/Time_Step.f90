!==============================================================================!
  subroutine Control_Mod_Time_Step(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real, intent(out) :: val
  logical, optional :: verbose
!==============================================================================!

  call Control % Read_Real_Item('TIME_STEP', 1.0e-2,      &
                                 val, verbose)

  end subroutine
