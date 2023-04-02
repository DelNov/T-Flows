!==============================================================================!
  subroutine Control_Mod_Simple_Underrelaxation_For_Momentum(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real, intent(out) :: val
  logical, optional :: verbose
!==============================================================================!

  call Control % Read_Real_Item('SIMPLE_UNDERRELAXATION_FOR_MOMENTUM', 0.6,  &
                                 val, verbose)

  end subroutine
