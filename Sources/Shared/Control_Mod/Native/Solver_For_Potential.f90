!==============================================================================!
  subroutine Control_Mod_Solver_For_Potential(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(SL), intent(out) :: val
  logical, optional          :: verbose
!==============================================================================!

  call Control_Mod_Read_Char_Item('SOLVER_FOR_POTENTIAL', 'cg',  &
                                   val, verbose)
  call String % To_Lower_Case(val)

  if( val .ne. 'bicg' .and. val .ne. 'cg') then
    call Message % Error(60,                                       &
             'Unknown linear solver for potential: '//trim(val)//  &
             '. \n This error is critical.  Exiting.',             &
             file=__FILE__, line=__LINE__, one_proc=.true.)
  end if

  end subroutine
