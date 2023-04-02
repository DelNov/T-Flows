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

  call Control % Read_Char_Item('TRACK_FRONT', 'no', val, verbose)
  call String % To_Upper_Case(val)

  if( val .eq. 'YES' ) then
    track_front = .true.

  else if( val .eq. 'NO' ) then
    track_front = .false.

  else
    call Message % Error(60,                                             &
                         'Unknown state for track front: '//trim(val)//  &
                         '. \n This error is critical.  Exiting.',       &
                         file=__FILE__, line=__LINE__, one_proc=.true.)
  end if

  end subroutine
