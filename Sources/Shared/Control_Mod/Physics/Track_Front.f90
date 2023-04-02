!==============================================================================!
  subroutine Track_Front(Control, track, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control
  logical, intent(out) :: track
  logical, optional    :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: val
!==============================================================================!

  call Control % Read_Char_Item('TRACK_FRONT', 'no', val, verbose)
  call String % To_Upper_Case(val)

  if( val .eq. 'YES' ) then
    track = .true.

  else if( val .eq. 'NO' ) then
    track = .false.

  else
    call Message % Error(60,                                             &
                         'Unknown state for track front: '//trim(val)//  &
                         '. \n This error is critical.  Exiting.',       &
                         file=__FILE__, line=__LINE__, one_proc=.true.)
  end if

  end subroutine
