!==============================================================================!
  subroutine Gu_Correction(Control, corr, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control
  logical, intent(out) :: corr
  logical, optional    :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: val
!==============================================================================!

  call Control % Read_Char_Item('GU_CORRECTION', 'no', val, verbose)
  call String % To_Upper_Case(val)

  if( val .eq. 'YES' ) then
    corr = .true.

  else if( val .eq. 'NO' ) then
    corr = .false.

  else
    call Message % Error(60,                                   &
             'Unknown state for Gu correction: '//trim(val)//  &
             '. \n This error is critical.  Exiting.',         &
             file=__FILE__, line=__LINE__, one_proc=.true.)
  end if

  end subroutine
