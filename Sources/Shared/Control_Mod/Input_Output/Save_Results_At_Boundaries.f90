!==============================================================================!
  subroutine Control_Mod_Save_Results_At_Boundaries(save_results_bnd, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical           :: save_results_bnd
  logical, optional :: verbose
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
