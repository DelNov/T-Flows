!==============================================================================!
  subroutine Control_Mod_Phase_Conductivities(val, verbose)
!------------------------------------------------------------------------------!
!   Reads as many densities as there are phases.                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real              :: val(0:1)
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  real :: def(2)
!==============================================================================!

  def = 1.0

  call Control_Mod_Read_Real_Array('PHASE_CONDUCTIVITIES', 2, def, val, verbose)

  end subroutine
