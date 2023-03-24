!==============================================================================!
  subroutine Control_Mod_Particle_Tracking(track_part, verbose)
!------------------------------------------------------------------------------!
!   Reading if Lagrangian particle tracking  will be used in simulations       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical, intent(out) :: track_part
  logical, optional    :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: val
!==============================================================================!

  call Control_Mod_Read_Char_Item('PARTICLE_TRACKING', 'no', val, verbose)
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
