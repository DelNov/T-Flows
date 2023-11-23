!==============================================================================!
  subroutine Phase_Capacities(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads phase thermal capacities for interface tracking VOF simulations.
!>  It is currently limited to two phases.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control   !! parent class
  real                :: val(0:1)  !! phase capacities of phases 0 and 1
  logical,   optional :: verbose   !! controls output verbosity
!-----------------------------------[Locals]-----------------------------------!
  real :: def(2)
!==============================================================================!

  def = 1.0

  call Control % Read_Real_Vector('PHASE_CAPACITIES', 2, def, val, verbose)

  end subroutine
