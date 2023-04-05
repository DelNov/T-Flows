!==============================================================================!
  subroutine Max_Least_Squares_Gradients_Iterations(Control, val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control
  integer, intent(out) :: val
  logical, optional    :: verbose
!==============================================================================!

  call Control % Read_Int_Item('MAX_LEAST_SQUARES_GRADIENTS_ITERATIONS', 4,  &
                                val, verbose)

  end subroutine
