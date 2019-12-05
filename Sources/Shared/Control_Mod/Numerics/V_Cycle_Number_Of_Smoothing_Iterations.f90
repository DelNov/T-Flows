!==============================================================================!
  subroutine Control_Mod_V_Cycle_Number_Of_Smoothing_Iterations(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer           :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Int_Item('V_CYCLE_NUMBER_OF_SMOOTHING_ITERATIONS',  &
                                  val, val, verbose)

  end subroutine
