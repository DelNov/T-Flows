!==============================================================================!
  subroutine Control_Mod_Distance_Function(d_func, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical           :: d_func
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(len=80) :: val
!==============================================================================!

  call Control_Mod_Read_Char_Item('DISTANCE_FUNCTION',  &
                                  'no', val, verbose)
  call To_Upper_Case(val)

  if( val .eq. 'YES' ) then
    d_func = .true.

  else if( val .eq. 'NO' ) then
    d_func = .false.

  else
    if(this_proc < 2) then
      print *, '# ERROR!  Unknown state for distance function: ', trim(val)
      print *, '# Exiting!'
    end if
    call Comm_Mod_End

  end if

  end subroutine
