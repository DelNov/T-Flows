!==============================================================================!
  subroutine Solver_For_Pressure(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads linear solver for pressure (or pressure correction) from control file.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)        :: Control  !! parent class
  character(SL), intent(out) :: val      !! linear solver (cg or bicg)
  logical, optional          :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Char_Item('SOLVER_FOR_PRESSURE', 'cg', val, verbose)
  call String % To_Lower_Case(val)

# if T_FLOWS_WITH_AMGX != 1
  if(val .eq. 'amgx') then
    call Message % Error(96,                                        &
             'You specified "amgx" solver for pressure but '    //  &
             'the code was compiled without the AMGX library! ' //  &
             '\n \n This error is critical.  Exiting.',             &
             file=__FILE__, line=__LINE__, one_proc=.true.)
  end if
# endif

  if(val .ne. 'bicg'      .and.  &
     val .ne. 'cg'        .and.  &
     val .ne. 'chebyshev' .and.  &
     val .ne. 'amgx') then
    call Message % Error(60,                                      &
             'Unknown linear solver for pressure: '//trim(val)//  &
             '. \n This error is critical.  Exiting.',            &
             file=__FILE__, line=__LINE__, one_proc=.true.)
  end if

  end subroutine
