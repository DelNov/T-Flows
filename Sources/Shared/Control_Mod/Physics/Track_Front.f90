!==============================================================================!
  subroutine Control_Mod_Track_Front(track_front, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical, intent(out) :: track_front
  logical, optional    :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: val
!==============================================================================!

  call Control_Mod_Read_Char_Item('TRACK_FRONT', 'no', val, verbose)
  call To_Upper_Case(val)

  if( val .eq. 'YES' ) then
    track_front = .true.

  else if( val .eq. 'NO' ) then
    track_front = .false.

  else
    if(this_proc < 2) then
      print *, '# ERROR!  Unknown state for track front: ',   &
                trim(val)
      print *, '# Exiting!'
    end if
    call Comm_Mod_End

  end if

  end subroutine
