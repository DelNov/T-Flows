!==============================================================================!
  subroutine Control_Mod_Phase_Densities(val, verbose)
!------------------------------------------------------------------------------!
!   Reads as many densities as there are phases.                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real, allocatable :: val(:)
  real              :: def(size(val))
  logical, optional :: verbose
!==============================================================================!

  def = 1.0

  call Control_Mod_Read_Real_Array('PHASE_DENSITIES',  &
                                    size(val), def, val, verbose)

  end subroutine
