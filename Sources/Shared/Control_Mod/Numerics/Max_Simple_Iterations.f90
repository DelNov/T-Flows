!==============================================================================!
  subroutine Max_Simple_Iterations(Control, val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control
  integer, intent(out) :: val
  logical, optional    :: verbose
!==============================================================================!

  call Control % Read_Int_Item('MAX_SIMPLE_ITERATIONS', 12,  &
                                val, verbose)

  end subroutine
