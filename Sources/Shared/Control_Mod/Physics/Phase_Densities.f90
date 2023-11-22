!==============================================================================!
  subroutine Phase_Densities(Control, val, verbose)
!------------------------------------------------------------------------------!
!   Reads as many densities as there are phases.                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control
  real                :: val(0:1)
  logical,   optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  real :: def(2)
!==============================================================================!

  def = 1.0

  call Control % Read_Real_Vector('PHASE_DENSITIES', 2, def, val, verbose)

  end subroutine
