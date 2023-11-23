!==============================================================================!
  subroutine Heat_Transfer(Control, heat, verbose)
!------------------------------------------------------------------------------!
!>  Reads, from the control file, if simulation involves heat transfer.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control  !! parent class
  logical, intent(out) :: heat     !! true if problem involves heat transfer
  logical, optional    :: verbose  !! controls output verbosity
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: val
!==============================================================================!

  call Control % Read_Char_Item('HEAT_TRANSFER', 'no', val, verbose)
  call String % To_Upper_Case(val)

  if( val .eq. 'YES' ) then
    heat = .true.

  else if( val .eq. 'NO' ) then
    heat = .false.

  else
    call Message % Error(60,                                               &
                         'Unknown state for heat transfer: '//trim(val)//  &
                         '. \n This error is critical.  Exiting.',         &
                         file=__FILE__, line=__LINE__, one_proc=.true.)
  end if

  end subroutine
