!==============================================================================!
  subroutine Control_Mod_Buoyancy(buoyancy, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical           :: buoyancy
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: val
!==============================================================================!

  call Control_Mod_Read_Char_Item('BUOYANCY', 'no', val, verbose)
  call To_Upper_Case(val)

  if( val .eq. 'YES' ) then
    buoyancy = .true.

  else if( val .eq. 'NO' ) then
    buoyancy = .false.

  else
    if(this_proc < 2) then
      print *, '# ERROR!  Unknown state for buoyancy: ', trim(val)
      print *, '# Exiting!'
    end if
    call Comm_Mod_End

  end if

  end subroutine
