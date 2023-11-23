!==============================================================================!
  subroutine Save_Results_At_Boundaries(Control, save_results_bnd, verbose)
!------------------------------------------------------------------------------!
!>  Reads, from the control file, if you should save results at boundarie.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control           !! parent class
  logical             :: save_results_bnd  !! should you save at boundaries
  logical,   optional :: verbose           !! controls output verbosity
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: val
!==============================================================================!

  call Control % Read_Char_Item('SAVE_RESULTS_AT_BOUNDARIES', 'yes',  &
                                 val, verbose)
  call String % To_Upper_Case(val)

  if( val .eq. 'YES' ) then
    save_results_bnd = .true.

  else if( val .eq. 'NO' ) then
    save_results_bnd = .false.

  else
    call Message % Error(72,                                                &
             'Unknown state for save results at boundaries: '//trim(val)//  &
             '. \n This error is critical.  Exiting.',                      &
             file=__FILE__, line=__LINE__, one_proc=.true.)
  end if

  end subroutine
