!==============================================================================!
  subroutine Control_Mod_Interface_Tracking(track_int, verbose)
!------------------------------------------------------------------------------!
!   Reading if vof will be used to model multiphase situation                  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical, intent(out) :: track_int
  logical, optional    :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: val
!==============================================================================!

  call Control_Mod_Read_Char_Item('INTERFACE_TRACKING', 'no', val, verbose)
  call To_Upper_Case(val)

  if( val .eq. 'YES' ) then
    track_int = .true.

  else if( val .eq. 'NO' ) then
    track_int = .false.

  else
    if(this_proc < 2) then
      print *, '# ERROR!  Unknown state for track front: ',   &
                trim(val)
      print *, '# Exiting!'
    end if
    call Comm_Mod_End

  end if

  end subroutine
