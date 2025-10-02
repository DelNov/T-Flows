!==============================================================================!
  subroutine Phase_Conductivities(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads phase thermal conductivites for interface tracking VOF simulations.
!>  It is currently limited to two phases.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control   !! parent class
  real                :: val(0:1)  !! phase conductivities of phases 0 and 1
  logical,   optional :: verbose   !! controls output verbosity
!-----------------------------------[Locals]-----------------------------------!
  real :: def(2)
!==============================================================================!

  def = 1.0

  call Control % Read_Real_Vector('PHASE_CONDUCTIVITIES', 2, def, val, verbose)

  end subroutine
