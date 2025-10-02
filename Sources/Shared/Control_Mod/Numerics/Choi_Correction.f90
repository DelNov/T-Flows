!==============================================================================!
  subroutine Choi_Correction(Control, corr, verbose)
!------------------------------------------------------------------------------!
!>  Reads from the control file if Choi correction will be used.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control  !! parent class
  logical, intent(out) :: corr     !! output value, true or false
  logical, optional    :: verbose  !! verbosity of the output
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: val
!==============================================================================!

  call Control % Read_Char_Item('CHOI_CORRECTION', 'no', val, verbose)
  call String % To_Upper_Case(val)

  if( val .eq. 'YES' ) then
    corr = .true.

  else if( val .eq. 'NO' ) then
    corr = .false.

  else
    call Message % Error(60,                                     &
             'Unknown state for Choi correction: '//trim(val)//  &
             '. \n This error is critical.  Exiting.',           &
             file=__FILE__, line=__LINE__, one_proc=.true.)
  end if

  end subroutine
