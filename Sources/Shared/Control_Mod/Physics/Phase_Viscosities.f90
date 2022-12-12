!==============================================================================!
  subroutine Control_Mod_Phase_Viscosities(val, verbose)
!------------------------------------------------------------------------------!
!   Reads as many viscosities as there are phases.                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real              :: val(0:1)
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  real :: def(2)
!==============================================================================!

  def = 1.0

  call Control_Mod_Read_Real_Array('PHASE_VISCOSITIES', 2, def, val, verbose)

  end subroutine
