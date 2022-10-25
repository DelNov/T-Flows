!==============================================================================!
  subroutine Control_Mod_Extrapolate_Temperature_Exp(temp_exp,  &
                                                     verbose)
!------------------------------------------------------------------------------!
!   Reading if temprature will be extrapolated to walls exponentially          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical, intent(out) :: temp_exp
  logical, optional    :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: val
!==============================================================================!

  call Control_Mod_Read_Char_Item('EXTRAPOLATE_TEMPERATURE_EXP',  &
                                  'no', val, verbose)
  call String % To_Upper_Case(val)

  if( val .eq. 'YES' ) then
    temp_exp = .true.

  else if( val .eq. 'NO' ) then
    temp_exp = .false.

  else
    if(this_proc < 2) then
      print *, '# ERROR!  Unknown state for temperature extrapolation: ',   &
                trim(val)
      print *, '# Exiting!'
    end if
    call Comm_Mod_End

  end if

  end subroutine
