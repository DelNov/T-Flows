!==============================================================================!
  subroutine Swarm_Coefficient_Of_Restitution(Control, cor, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control
  real,   intent(out) :: cor
  logical,   optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  real :: def
!==============================================================================!

  def = 1.0

  call Control % Read_Real_Item('SWARM_COEFFICIENT_OF_RESTITUTION',  &
                                 def, cor, verbose)

  end subroutine
