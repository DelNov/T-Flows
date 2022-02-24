!==============================================================================!
  subroutine Control_Mod_Phase_Viscosities(val, verbose)
!------------------------------------------------------------------------------!
!   Reads as many viscosities as there are phases.                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real, allocatable :: val(:)
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  real :: def(size(val))
!==============================================================================!

  def = 1.0

  call Control_Mod_Read_Real_Array('PHASE_VISCOSITIES',  &
                                    size(val), def, val, verbose)

  end subroutine
