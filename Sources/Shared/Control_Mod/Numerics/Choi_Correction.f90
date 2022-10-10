!==============================================================================!
  subroutine Control_Mod_Choi_Correction(choi_correction, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical, intent(out) :: choi_correction
  logical, optional    :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: val
!==============================================================================!

  call Control_Mod_Read_Char_Item('CHOI_CORRECTION', 'no', val, verbose)
  call String % To_Upper_Case(val)

  if( val .eq. 'YES' ) then
    choi_correction = .true.

  else if( val .eq. 'NO' ) then
    choi_correction = .false.

  else
    if(this_proc < 2) then
      print *, '# ERROR!  Unknown state for Choi''s correction: ',   &
                trim(val)
      print *, '# Exiting!'
    end if
    call Comm_Mod_End
    stop

  end if

  end subroutine
