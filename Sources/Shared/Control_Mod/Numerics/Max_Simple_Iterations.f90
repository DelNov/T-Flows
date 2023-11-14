!==============================================================================!
  subroutine Max_Simple_Iterations(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads maximum number of SIMPLE iterations.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control  !! parent class
  integer, intent(out) :: val      !! max SIMPLE iterations
  logical, optional    :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Int_Item('MAX_SIMPLE_ITERATIONS', 12,  &
                                val, verbose)

  end subroutine
