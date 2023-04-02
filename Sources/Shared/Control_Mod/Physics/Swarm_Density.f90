!==============================================================================!
  subroutine Swarm_Density(Control, s_den, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control
  real,   intent(out) :: s_den
  logical,   optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  real :: def
!==============================================================================!

  def = 1000.0

  call Control % Read_Real_Item('SWARM_DENSITY', def, s_den, verbose)

  end subroutine
