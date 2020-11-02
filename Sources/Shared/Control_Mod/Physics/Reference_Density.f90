!==============================================================================!
  subroutine Control_Mod_Reference_Density(d_ref, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real              :: d_ref
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  real :: def
!==============================================================================!

  def = 1000.0

  call Control_Mod_Read_Real_Item('REFERENCE_DENSITY', def, d_ref, verbose)

  end subroutine