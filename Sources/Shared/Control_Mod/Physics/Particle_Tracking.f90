!==============================================================================!
  subroutine Particle_Tracking(Control, track_part, verbose)
!------------------------------------------------------------------------------!
!>  Reads from control file if Lagrangian particle tracking  will be used in
!>  simulations.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control     !! parent class
  logical, intent(out) :: track_part  !! true to engage particle tracking
  logical, optional    :: verbose     !! controls output verbosity
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: val
!==============================================================================!

  call Control % Read_Char_Item('PARTICLE_TRACKING', 'no', val, verbose)
  call String % To_Upper_Case(val)

  if( val .eq. 'YES' ) then
    track_part = .true.

  else if( val .eq. 'NO' ) then
    track_part = .false.

  else
    call Message % Error(72,                                                &
                      'Unknown state for particle tracking: '//trim(val)//  &
                      '. \n This error is critical.  Exiting.',             &
                      file=__FILE__, line=__LINE__, one_proc=.true.)
  end if

  end subroutine
