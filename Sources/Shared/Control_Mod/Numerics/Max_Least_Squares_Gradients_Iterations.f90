!==============================================================================!
  subroutine Max_Least_Squares_Gradients_Iterations(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads maximum number of iterations for least-squares gradient method.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control  !! parent class
  integer, intent(out) :: val      !! max least-squares iteration
  logical, optional    :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Int_Item('MAX_LEAST_SQUARES_GRADIENTS_ITERATIONS', 4,  &
                                val, verbose)

  end subroutine
