!==============================================================================!
  subroutine Control_Mod_Save_Results_Units(save_results_units, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical           :: save_results_units
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: val
!==============================================================================!

  call Control_Mod_Read_Char_Item('SAVE_RESULTS_UNITS', 'yes',  &
                                   val, verbose)
  call String % To_Upper_Case(val)

  if( val .eq. 'YES' ) then
    save_results_units = .true.

  else if( val .eq. 'NO' ) then
    save_results_units = .false.

  else
    if(this_proc < 2) then
      print *, '# ERROR!  Unknown state for save_results_units: ', trim(val)
      print *, '# Exiting!'
    end if
    call Comm_Mod_End

  end if

  end subroutine
