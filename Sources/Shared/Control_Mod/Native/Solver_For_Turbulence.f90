!==============================================================================!
  subroutine Solver_For_Turbulence(Control, val, verbose)
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

  call Control % Read_Char_Item('SOLVER_FOR_TURBULENCE', 'bicg', val, verbose)
  call String % To_Lower_Case(val)

  if(val .ne. 'bicg' .and. val .ne. 'cg') then
    call Message % Error(72,                                             &
             'Unknown linear solver for wall turbulence: '//trim(val)//  &
             '. \n This error is critical.  Exiting.',                   &
             file=__FILE__, line=__LINE__, one_proc=.true.)
  end if

  end subroutine
