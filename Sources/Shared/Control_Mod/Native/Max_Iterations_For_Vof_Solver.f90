!==============================================================================!
  subroutine Max_Iterations_For_Vof_Solver(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads maximum number of iterations (in a linear solver) for VOF function
!>  from the control file.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control  !! parent class
  integer, intent(out) :: val      !! max iterations
  logical, optional    :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Int_Item('MAX_ITERATIONS_FOR_VOF_SOLVER',  &
                                6, val, verbose)

  end subroutine
