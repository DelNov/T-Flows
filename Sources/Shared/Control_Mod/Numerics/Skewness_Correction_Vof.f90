!==============================================================================!
  subroutine Control_Mod_Skewness_Correction_Vof(skew_corr, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical, intent(out) :: skew_corr
  logical, optional    :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: val
!==============================================================================!

  call Control_Mod_Read_Char_Item('SKEWNESS_CORRECTION_VOF',   &
                                  'no', val, verbose)
  call String % To_Upper_Case(val)

  if( val .eq. 'YES' ) then
    skew_corr = .true.

  else if( val .eq. 'NO' ) then
    skew_corr = .false.

  else
    call Message % Error(72,                                         &
             'Unknown state for skewness correction: '//trim(val)//  &
             '. \n This error is critical.  Exiting.',               &
             file=__FILE__, line=__LINE__, one_proc=.true.)
  end if

  end subroutine
