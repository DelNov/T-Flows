!==============================================================================!
  subroutine Min_Simple_Iterations(Control, val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control
  integer, intent(out) :: val
  logical, optional    :: verbose
!==============================================================================!

  call Control % Read_Int_Item('MIN_SIMPLE_ITERATIONS', 3,  &
                                val, verbose)

  end subroutine
