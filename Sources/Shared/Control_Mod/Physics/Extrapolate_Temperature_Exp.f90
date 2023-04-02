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

  call Control % Read_Char_Item('EXTRAPOLATE_TEMPERATURE_EXP',  &
                                'no', val, verbose)
  call String % To_Upper_Case(val)

  if( val .eq. 'YES' ) then
    temp_exp = .true.

  else if( val .eq. 'NO' ) then
    temp_exp = .false.

  else
    call Message % Error(72,                                               &
             'Unknown state for temperature extrapolation: '//trim(val)//  &
             '. \n This error is critical.  Exiting.',                     &
             file=__FILE__, line=__LINE__, one_proc=.true.)
  end if

  end subroutine
