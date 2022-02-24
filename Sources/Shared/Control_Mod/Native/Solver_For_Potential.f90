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
  call To_Lower_Case(val)

  if( val.ne.'bicg' .and. val.ne.'cgs' .and. val.ne.'cg') then
    if(this_proc < 2) then
      print *, '# ERROR!  Unknown linear solver for pressure: ', trim(val)
      print *, '# Exiting!'
    end if
    call Comm_Mod_End
  end if

  end subroutine
