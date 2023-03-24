!==============================================================================!
  subroutine Control_Mod_Gu_Correction(gu_correction, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical, intent(out) :: gu_correction
  logical, optional    :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: val
!==============================================================================!

  call Control_Mod_Read_Char_Item('GU_CORRECTION', 'no', val, verbose)
  call String % To_Upper_Case(val)

  if( val .eq. 'YES' ) then
    gu_correction = .true.

  else if( val .eq. 'NO' ) then
    gu_correction = .false.

  else
    call Message % Error(60,                                   &
             'Unknown state for Gu correction: '//trim(val)//  &
             '. \n This error is critical.  Exiting.',         &
             file=__FILE__, line=__LINE__, one_proc=.true.)
  end if

  end subroutine
