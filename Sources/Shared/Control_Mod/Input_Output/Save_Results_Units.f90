!==============================================================================!
  subroutine Save_Results_Units(Control, save_results_uni, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control
  logical             :: save_results_uni
  logical, optional   :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: val
!==============================================================================!

  call Control % Read_Char_Item('SAVE_RESULTS_UNITS', 'yes',  &
                                 val, verbose)
  call String % To_Upper_Case(val)

  if( val .eq. 'YES' ) then
    save_results_uni = .true.

  else if( val .eq. 'NO' ) then
    save_results_uni = .false.

  else
    call Message % Error(72,                                                &
             'Unknown state for save results units: '//trim(val)//          &
             '. \n This error is critical.  Exiting.',                      &
             file=__FILE__, line=__LINE__, one_proc=.true.)
  end if

  end subroutine
