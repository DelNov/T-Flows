!==============================================================================!
  subroutine Control_Mod_Phase_Change(phase_change, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical           :: phase_change
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: val
!==============================================================================!

  call Control_Mod_Read_Char_Item('PHASE_CHANGE', 'no', val, verbose)
  call To_Upper_Case(val)

  if( val .eq. 'YES' ) then
    phase_change = .true.

  else if( val .eq. 'NO' ) then
    phase_change = .false.

  else
    if(this_proc < 2) then
      print *, '# ERROR!  Unknown state for phase_change: ', trim(val)
      print *, '# Exiting!'
    end if
    call Comm_Mod_End

  end if

  end subroutine
