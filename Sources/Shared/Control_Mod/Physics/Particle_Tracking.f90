!==============================================================================!
  subroutine Control_Mod_Particle_Tracking(track_part, verbose)
!------------------------------------------------------------------------------!
!   Reading if Lagrangian particle tracking  will be used in simulations       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical           :: track_part
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: val
!==============================================================================!

  call Control_Mod_Read_Char_Item('PARTICLE_TRACKING', 'no', val, verbose)
  call To_Upper_Case(val)

  if( val .eq. 'YES' ) then
    track_part = .true.

  else if( val .eq. 'NO' ) then
    track_part = .false.

  else
    if(this_proc < 2) then
      print *, '# ERROR!  Unknown state for track front: ',   &
                trim(val)
      print *, '# Exiting!'
    end if
    call Comm_Mod_End

  end if

  end subroutine
