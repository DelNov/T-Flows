!==============================================================================!
  subroutine Track_Front(Control, track, verbose)
!------------------------------------------------------------------------------!
!>  Reads if front will be re-contstructed with VOF simulations of a multiphase
!>  flow. Front is represented with polygons spanned inside computational cells.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control  !! parent class
  logical, intent(out) :: track    !! true if front is tracked
  logical, optional    :: verbose  !! controls output verbosity
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
