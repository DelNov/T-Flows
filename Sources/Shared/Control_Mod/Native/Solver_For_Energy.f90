!==============================================================================!
  subroutine Control_Mod_Solver_For_Energy(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(SL), intent(out) :: val
  logical, optional          :: verbose
!==============================================================================!

  call Control % Read_Char_Item('SOLVER_FOR_ENERGY', 'cg', val, verbose)
  call String % To_Lower_Case(val)

  if( val .ne. 'bicg' .and. val .ne. 'cg') then
    call Message % Error(60,                                    &
             'Unknown linear solver for energy: '//trim(val)//  &
             '. \n This error is critical.  Exiting.',          &
             file=__FILE__, line=__LINE__, one_proc=.true.)
  end if

  end subroutine
