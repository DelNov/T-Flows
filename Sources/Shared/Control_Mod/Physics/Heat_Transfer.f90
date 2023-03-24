!==============================================================================!
  subroutine Control_Mod_Heat_Transfer(heat_transfer, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical, intent(out) :: heat_transfer
  logical, optional    :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: val
!==============================================================================!

  call Control_Mod_Read_Char_Item('HEAT_TRANSFER', 'no', val, verbose)
  call String % To_Upper_Case(val)

  if( val .eq. 'YES' ) then
    heat_transfer = .true.

  else if( val .eq. 'NO' ) then
    heat_transfer = .false.

  else
    call Message % Error(60,                                               &
                         'Unknown state for heat transfer: '//trim(val)//  &
                         '. \n This error is critical.  Exiting.',         &
                         file=__FILE__, line=__LINE__, one_proc=.true.)
  end if

  end subroutine
