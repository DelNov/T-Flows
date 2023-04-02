!==============================================================================!
  subroutine Control_Mod_Choi_Correction(choi_correction, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical, intent(out) :: choi_correction
  logical, optional    :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: val
!==============================================================================!

  call Control % Read_Char_Item('CHOI_CORRECTION', 'no', val, verbose)
  call String % To_Upper_Case(val)

  if( val .eq. 'YES' ) then
    choi_correction = .true.

  else if( val .eq. 'NO' ) then
    choi_correction = .false.

  else
    call Message % Error(60,                                     &
             'Unknown state for Choi correction: '//trim(val)//  &
             '. \n This error is critical.  Exiting.',           &
             file=__FILE__, line=__LINE__, one_proc=.true.)
  end if

  end subroutine
