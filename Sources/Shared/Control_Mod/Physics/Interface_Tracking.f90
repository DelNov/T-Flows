!==============================================================================!
  subroutine Interface_Tracking(Control, track_int, verbose)
!------------------------------------------------------------------------------!
!   Reading if vof will be used to model multiphase situation                  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control
  logical, intent(out) :: track_int
  logical, optional    :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: val
!==============================================================================!

  call Control % Read_Char_Item('INTERFACE_TRACKING', 'no', val, verbose)
  call String % To_Upper_Case(val)

  if( val .eq. 'YES' ) then
    track_int = .true.

  else if( val .eq. 'NO' ) then
    track_int = .false.

  else
    call Message % Error(72,                                                 &
                      'Unknown state for interface tracking: '//trim(val)//  &
                      '. \n This error is critical.  Exiting.',              &
                      file=__FILE__, line=__LINE__, one_proc=.true.)
  end if

  end subroutine
