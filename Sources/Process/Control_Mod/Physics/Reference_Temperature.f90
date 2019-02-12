!==============================================================================!
  subroutine Control_Mod_Reference_Temperature(t_ref, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real              :: t_ref
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  real :: def
!==============================================================================!

  def = 20.0

  call Control_Mod_Read_Real_Item('REFERENCE_TEMPERATURE', def, t_ref, verbose)

  end subroutine
