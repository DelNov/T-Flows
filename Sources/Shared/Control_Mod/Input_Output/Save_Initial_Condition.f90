!==============================================================================!
  subroutine Control_Mod_Save_Initial_Condition(save_init_cond, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical           :: save_init_cond
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: val
!==============================================================================!

  call Control_Mod_Read_Char_Item('SAVE_INITIAL_CONDITION', 'yes',  &
                                   val, verbose)
  call To_Upper_Case(val)

  if( val .eq. 'YES' ) then
    save_init_cond = .true.

  else if( val .eq. 'NO' ) then
    save_init_cond = .false.

  else
    if(this_proc < 2) then
      print *, '# ERROR!  Unknown state for save_init_cond: ', trim(val)
      print *, '# Exiting!'
    end if
    call Comm_Mod_End

  end if

  end subroutine
