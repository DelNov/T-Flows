!==============================================================================!
  subroutine Mass_Transfer(Control, phase_change, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control
  logical, intent(out) :: phase_change
  logical, optional    :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: val
!==============================================================================!

  call Control % Read_Char_Item('MASS_TRANSFER', 'no', val, verbose)
  call String % To_Upper_Case(val)

  if( val .eq. 'YES' ) then
    phase_change = .true.

  else if( val .eq. 'NO' ) then
    phase_change = .false.

  else
    call Message % Error(60,                                              &
                         'Unknown state for phase change: '//trim(val)//  &
                         '. \n This error is critical.  Exiting.',        &
                         file=__FILE__, line=__LINE__, one_proc=.true.)
  end if

  end subroutine
