!==============================================================================!
  subroutine Control_Mod_Simple_Underrelaxation_For_Energy(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real, intent(out) :: val
  logical, optional :: verbose
!==============================================================================!

  call Control % Read_Real_Item('SIMPLE_UNDERRELAXATION_FOR_ENERGY', 0.5,  &
                                 val, verbose)

  end subroutine
