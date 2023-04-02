!==============================================================================!
  subroutine Control_Mod_Swarm_Diameter(s_dia, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real, intent(out) :: s_dia
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  real :: def
!==============================================================================!

  def = 2.5e-5

  call Control % Read_Real_Item('SWARM_DIAMETER', def, s_dia, verbose)

  end subroutine
