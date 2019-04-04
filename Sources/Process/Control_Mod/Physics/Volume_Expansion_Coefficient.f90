!==============================================================================!
  subroutine Control_Mod_Volume_Expansion_Coefficient(cor, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real              :: cor
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  real :: def
!==============================================================================!

  def = 1.0

  call Control_Mod_Read_Real_Item('VOLUME_EXPANSION_COEFFICIENT',  &
                                   def, cor, verbose)

  end subroutine
