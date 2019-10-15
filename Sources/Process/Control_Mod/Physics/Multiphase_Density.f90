!==============================================================================!
  subroutine Control_Mod_Multiphase_Density(val, verbose)
!------------------------------------------------------------------------------!
!   Reading stuff related to passive scalars                                   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real, allocatable :: val(:)
  real              :: def(size(val))
  logical, optional :: verbose
!==============================================================================!

  def = 1.0

  call Control_Mod_Read_Real_Array('DENSITY_MULT', size(val), def, val, verbose)

  end subroutine
