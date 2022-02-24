!==============================================================================!
  subroutine Control_Mod_Max_Least_Squares_Gradients_Iterations(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer, intent(out) :: val
  logical, optional    :: verbose
!==============================================================================!

  call Control_Mod_Read_Int_Item('MAX_LEAST_SQUARES_GRADIENTS_ITERATIONS', 4,  &
                                  val, verbose)

  end subroutine
