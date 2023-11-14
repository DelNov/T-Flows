!==============================================================================!
  subroutine Report_Volume_Balance(Control, vol_bal, verbose)
!------------------------------------------------------------------------------!
!>  Reads the control file to find out if volume balance should be reported.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control  !! parent class
  logical, intent(out) :: vol_bal  !! report balance (true or false)
  logical, optional    :: verbose  !! controls output verbosity
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: val
!==============================================================================!

  call Control % Read_Char_Item('REPORT_VOLUME_BALANCE',   &
                                'no', val, verbose)
  call String % To_Upper_Case(val)

  if( val .eq. 'YES' ) then
    vol_bal = .true.

  else if( val .eq. 'NO' ) then
    vol_bal = .false.

  else
    call Message % Error(72,                                              &
             'Unknown state for reporting volume balance: '//trim(val)//  &
             '. \n This error is critical.  Exiting.',                    &
             file=__FILE__, line=__LINE__, one_proc=.true.)
  end if

  end subroutine
