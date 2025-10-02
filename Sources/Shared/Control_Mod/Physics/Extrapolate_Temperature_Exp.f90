!==============================================================================!
  subroutine Extrapolate_Temperature_Exp(Control, temp_exp, verbose)
!------------------------------------------------------------------------------!
!>  Reads, from the control file, if temprature will be extrapolated to the
!>  walls with exponential interpolation.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control   !! parent class
  logical, intent(out) :: temp_exp  !! true if temperatures are extrapolated
                                    !! exponentially, false otherwise
  logical, optional    :: verbose   !! controls output verbosity
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
