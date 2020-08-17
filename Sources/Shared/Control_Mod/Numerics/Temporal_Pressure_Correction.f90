!==============================================================================!
  subroutine Control_Mod_Temporal_Pressure_Correction(temp_corr, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical           :: temp_corr
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: val
!==============================================================================!

  call Control_Mod_Read_Char_Item('TEMPORAL_PRESSURE_CORRECTION',   &
                                  'no', val, verbose)
  call To_Upper_Case(val)

  if( val .eq. 'YES' ) then
    temp_corr = .true.

  else if( val .eq. 'NO' ) then
    temp_corr = .false.

  else
    if(this_proc < 2) then
      print *, '# ERROR!  Unknown state for temporal pressure correction: ',   &
                trim(val)
      print *, '# Exiting!'
    end if
    call Comm_Mod_End

  end if

  end subroutine
