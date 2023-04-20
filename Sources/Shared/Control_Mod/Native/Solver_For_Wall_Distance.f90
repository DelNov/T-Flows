!==============================================================================!
  subroutine Solver_For_Wall_Distance(Control, val, verbose)
!------------------------------------------------------------------------------!
!   Default values in this, and her sister functions, is simply as follows:    !
!   If matrix is expected to be symmetric, take 'cg', otherwise take 'bicg'.   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)        :: Control
  character(SL), intent(out) :: val
  logical, optional          :: verbose
!==============================================================================!

  call Control % Read_Char_Item('SOLVER_FOR_WALL_DISTANCE', 'cg', val, verbose)
  call String % To_Lower_Case(val)

  if(val .ne. 'bicg' .and. val .ne. 'cg') then
    call Message % Error(60,                                           &
             'Unknown linear solver for wall distance: '//trim(val)//  &
             '. \n This error is critical.  Exiting.',                 &
             file=__FILE__, line=__LINE__, one_proc=.true.)
  end if

  end subroutine
