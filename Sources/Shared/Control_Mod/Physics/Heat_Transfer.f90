!==============================================================================!
  subroutine Control_Mod_Heat_Transfer(heat_transfer, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical           :: heat_transfer
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(len=80) :: val
!==============================================================================!

  call Control_Mod_Read_Char_Item('HEAT_TRANSFER', 'no', val, verbose)
  call To_Upper_Case(val)

  if( val .eq. 'YES' ) then
    heat_transfer = .true.

  else if( val .eq. 'NO' ) then
    heat_transfer = .false.

  else
    if(this_proc < 2) then
      print *, '# ERROR!  Unknown state for heat transfer: ', trim(val)
      print *, '# Exiting!'
    end if
    call Comm_Mod_End

  end if

  end subroutine
