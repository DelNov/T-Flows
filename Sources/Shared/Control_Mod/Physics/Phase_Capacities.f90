!==============================================================================!
  subroutine Control_Mod_Phase_Capacities(val, verbose)
!------------------------------------------------------------------------------!
!   Reads as many capacities as there are phases.                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real, allocatable :: val(:)
  real              :: def(size(val))
  logical, optional :: verbose
!==============================================================================!

  def = 1.0

  call Control_Mod_Read_Real_Array('PHASE_CAPACITIES',  &
                                    size(val), def, val, verbose)

  end subroutine
