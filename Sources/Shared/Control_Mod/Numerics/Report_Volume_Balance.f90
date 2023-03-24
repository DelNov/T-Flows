!==============================================================================!
  subroutine Control_Mod_Report_Volume_Balance(vol_bal, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical, intent(out) :: vol_bal
  logical, optional    :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: val
!==============================================================================!

  call Control_Mod_Read_Char_Item('REPORT_VOLUME_BALANCE',   &
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
