!==============================================================================!
  subroutine Control_Mod_Phase_Capacities(val, verbose)
!------------------------------------------------------------------------------!
!   Reads as many capacities as there are phases.                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real              :: val(0:1)
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  real :: def(2)
!==============================================================================!

  def = 1.0

  call Control % Read_Real_Vector('PHASE_CAPACITIES', 2, def, val, verbose)

  end subroutine
