!==============================================================================!
  subroutine Mass_Transfer(Control, phase_change, verbose)
!------------------------------------------------------------------------------!
!>  Reads if the simulation involves mass transfer (phase change) from the
!>  control file.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control       !! parent class
  logical, intent(out) :: phase_change  !! true if phase change is modelled
  logical, optional    :: verbose       !! controls output verbosity
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
