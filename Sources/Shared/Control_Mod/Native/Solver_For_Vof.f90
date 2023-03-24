!==============================================================================!
  subroutine Control_Mod_Solver_For_Vof(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(SL), intent(out) :: val
  logical, optional          :: verbose
!==============================================================================!

  call Control_Mod_Read_Char_Item('SOLVER_FOR_VOF', 'cg',  &
                                   val, verbose)
  call String % To_Lower_Case(val)

  if( val .ne. 'bicg' .and. val .ne. 'cg') then
    call Message % Error(72,                                             &
             'Unknown linear solver for wall multiphase: '//trim(val)//  &
             '. \n This error is critical.  Exiting.',                   &
             file=__FILE__, line=__LINE__, one_proc=.true.)
  end if

  end subroutine
