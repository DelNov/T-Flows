!==============================================================================!
  subroutine Max_Gauss_Gradients_Iterations(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads maximum number of iterations for Gaussian gradient method.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control  !! parent class
  integer, intent(out) :: val      !! max Gauss-gradient iterations
  logical, optional    :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Int_Item('MAX_GAUSS_GRADIENTS_ITERATIONS', 12,  &
                                val, verbose)

  end subroutine
