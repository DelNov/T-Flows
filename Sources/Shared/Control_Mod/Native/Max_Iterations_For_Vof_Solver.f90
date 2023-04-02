!==============================================================================!
  subroutine Max_Iterations_For_Vof_Solver(Control, val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control
  integer, intent(out) :: val
  logical, optional    :: verbose
!==============================================================================!

  call Control % Read_Int_Item('MAX_ITERATIONS_FOR_VOF_SOLVER',  &
                                6, val, verbose)

  end subroutine
