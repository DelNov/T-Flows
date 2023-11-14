!==============================================================================!
  subroutine Track_Surface(Control, track_surf, verbose)
!------------------------------------------------------------------------------!
!>  Reads if surface will be re-contstructed with VOF simulations of a
!>  multiphase flow. Surface is represented with triangular grids independent
!>  from computational cells.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control     !! parent class
  logical, intent(out) :: track_surf  !! true if surface is tracked
  logical, optional    :: verbose     !! controls output verbosity
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
