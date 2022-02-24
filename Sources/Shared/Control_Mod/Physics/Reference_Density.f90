!==============================================================================!
  subroutine Control_Mod_Reference_Density(d_ref, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real, intent(out) :: d_ref
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  real :: def
!==============================================================================!

  def = 0.0

  call Control_Mod_Read_Real_Item('REFERENCE_DENSITY', def, d_ref, verbose)

  end subroutine
