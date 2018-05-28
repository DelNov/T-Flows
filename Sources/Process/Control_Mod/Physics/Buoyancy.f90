!==============================================================================!
  subroutine Control_Mod_Buoyancy(verbose)
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Flow_Mod,  only: buoyancy
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(len=80) :: val
!==============================================================================!

  call Control_Mod_Read_Char_Item('BUOYANCY', 'no',  &
                                   val, verbose)
  call To_Upper_Case(val)

  if( val .eq. 'YES' ) then
    buoyancy = YES

  else if( val .eq. 'NO' ) then
    buoyancy = NO

  else
    print *, '# Unknown state for buoyancy: ', trim(val)
    print *, '# Exiting!'
    stop

  end if

  end subroutine
