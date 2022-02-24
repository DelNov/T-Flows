!==============================================================================!
  subroutine Control_Mod_Swarm_Coefficient_Of_Restitution(cor, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real, intent(out) :: cor
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  real :: def
!==============================================================================!

  def = 1.0

  call Control_Mod_Read_Real_Item('SWARM_COEFFICIENT_OF_RESTITUTION',  &
                                   def, cor, verbose)

  end subroutine
