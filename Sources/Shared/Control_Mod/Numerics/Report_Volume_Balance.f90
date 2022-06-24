!==============================================================================!
  subroutine Control_Mod_Report_Volume_Balance(vol_bal, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical           :: vol_bal
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: val
!==============================================================================!

  call Control_Mod_Read_Char_Item('REPORT_VOLUME_BALANCE',   &
                                  'no', val, verbose)
  call To_Upper_Case(val)

  if( val .eq. 'YES' ) then
    vol_bal = .true.

  else if( val .eq. 'NO' ) then
    vol_bal = .false.

  else
    if(this_proc < 2) then
      print *, '# ERROR!  Unknown state for reporting volume balance: ',   &
                trim(val)
      print *, '# Exiting!'
    end if
    call Comm_Mod_End

  end if

  end subroutine
