!==============================================================================!
  subroutine Max_Gauss_Gradients_Iterations(Control, val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control
  integer, intent(out) :: val
  logical, optional    :: verbose
!==============================================================================!

  call Control % Read_Int_Item('MAX_GAUSS_GRADIENTS_ITERATIONS', 12,  &
                                val, verbose)

  end subroutine
