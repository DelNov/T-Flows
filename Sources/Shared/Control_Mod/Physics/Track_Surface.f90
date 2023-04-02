!==============================================================================!
  subroutine Control_Mod_Track_Surface(track_surf, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical, intent(out) :: track_surf
  logical, optional    :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: val
!==============================================================================!

  call Control % Read_Char_Item('TRACK_SURFACE', 'no', val, verbose)
  call String % To_Upper_Case(val)

  if( val .eq. 'YES' ) then
    track_surf = .true.

  else if( val .eq. 'NO' ) then
    track_surf = .false.

  else
    call Message % Error(60,                                               &
                         'Unknown state for track surface: '//trim(val)//  &
                         '. \n This error is critical.  Exiting.',         &
                         file=__FILE__, line=__LINE__, one_proc=.true.)
  end if

  end subroutine
