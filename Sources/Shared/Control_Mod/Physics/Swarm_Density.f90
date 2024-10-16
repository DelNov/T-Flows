!==============================================================================!
  subroutine Control_Mod_Swarm_Density(s_den, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real, intent(out) :: s_den
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  real :: def
!==============================================================================!

  def = 1000.0

  call Control_Mod_Read_Real_Item('SWARM_DENSITY', def, s_den, verbose)

  end subroutine
