!==============================================================================!
  subroutine Solver_For_Energy(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads linear solver for energy (enthalpy) from control file.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)        :: Control  !! parent class
  character(SL), intent(out) :: val      !! linear solver (cg or bicg)
  logical, optional          :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Char_Item('SOLVER_FOR_ENERGY', 'bicg', val, verbose)
  call String % To_Lower_Case(val)

  if(val .ne. 'bicg' .and. val .ne. 'cg') then
    call Message % Error(60,                                    &
             'Unknown linear solver for energy: '//trim(val)//  &
             '. \n This error is critical.  Exiting.',          &
             file=__FILE__, line=__LINE__, one_proc=.true.)
  end if

  end subroutine
