!==============================================================================!
  subroutine Control_Mod_V_Cycle_Residual_Ratio(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real              :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Real_Item('V_CYCLE_RESIDUAL_RATIO', &
                                   val, val, verbose)

  end subroutine
