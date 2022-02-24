!==============================================================================!
  subroutine Control_Mod_Phase_Conductivities(val, verbose)
!------------------------------------------------------------------------------!
!   Reads as many densities as there are phases.                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real, allocatable :: val(:)
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  real :: def(size(val))
!==============================================================================!

  def = 1.0

  call Control_Mod_Read_Real_Array('PHASE_CONDUCTIVITIES',  &
                                    size(val), def, val, verbose)

  end subroutine
