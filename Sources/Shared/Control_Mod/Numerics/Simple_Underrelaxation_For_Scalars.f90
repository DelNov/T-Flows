!==============================================================================!
  subroutine Control_Mod_Simple_Underrelaxation_For_Scalars(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real, intent(out) :: val
  logical, optional :: verbose
!==============================================================================!

  call Control % Read_Real_Item('SIMPLE_UNDERRELAXATION_FOR_SCALARS',  &
                                 0.5, val, verbose)

  end subroutine
