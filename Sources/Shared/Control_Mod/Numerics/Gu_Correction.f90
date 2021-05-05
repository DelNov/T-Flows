!==============================================================================!
  subroutine Control_Mod_Gu_Correction(gu_correction, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical           :: gu_correction
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: val
!==============================================================================!

  call Control_Mod_Read_Char_Item('GU_CORRECTION', 'yes', val, verbose)
  call To_Upper_Case(val)

  if( val .eq. 'YES' ) then
    gu_correction = .true.

  else if( val .eq. 'NO' ) then
    gu_correction = .false.

  else
    if(this_proc < 2) then
      print *, '# ERROR!  Unknown state for track front: ',   &
                trim(val)
      print *, '# Exiting!'
    end if
    call Comm_Mod_End
    stop

  end if

  end subroutine
