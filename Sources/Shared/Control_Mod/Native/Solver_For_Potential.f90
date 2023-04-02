!==============================================================================!
  subroutine Solver_For_Potential(Control, val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)        :: Control
  character(SL), intent(out) :: val
  logical, optional          :: verbose
!==============================================================================!

  call Control % Read_Char_Item('SOLVER_FOR_POTENTIAL', 'cg', val, verbose)
  call String % To_Lower_Case(val)

  if( val .ne. 'bicg' .and. val .ne. 'cg') then
    call Message % Error(60,                                       &
             'Unknown linear solver for potential: '//trim(val)//  &
             '. \n This error is critical.  Exiting.',             &
             file=__FILE__, line=__LINE__, one_proc=.true.)
  end if

  end subroutine
