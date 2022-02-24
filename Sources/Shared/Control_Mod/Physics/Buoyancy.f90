!==============================================================================!
  subroutine Control_Mod_Buoyancy(val, verbose)
!------------------------------------------------------------------------------!
!   Reading buoyancy model from the control file.                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(SL), intent(out) :: val
  logical,       optional    :: verbose
!==============================================================================!

  call Control_Mod_Read_Char_Item('BUOYANCY', 'none', val, verbose)

  call To_Upper_Case(val)

  end subroutine
