!==============================================================================!
  subroutine Control_Mod_Solver_For_Pressure(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(len=80) :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Char_Item('SOLVER_FOR_PRESSURE', 'bicg',  &
                                   val, verbose)
  call To_Upper_Case(val)

  if( val.ne.'BICG' .and. val.ne.'CGS'  .and. val.ne.'CG') then
    print *, '# Unknown linear solver for pressure: ', trim(val)
    print *, '# Exiting!'
    stop 
  end if

  end subroutine
