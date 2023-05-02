!==============================================================================!
  subroutine Swarm_Diameter(Control, s_dia, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control
  real,   intent(out) :: s_dia
  logical,   optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  real :: def
!==============================================================================!

  def = 2.5e-5

  call Control % Read_Real_Item('SWARM_DIAMETER', def, s_dia, verbose)

  end subroutine
