!==============================================================================!
  subroutine Load_Physical_Properties(grid)
!------------------------------------------------------------------------------!
!   Reads physical properties from control file.                               !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
  use Flow_Mod
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!==============================================================================!

  call Control_Mod_Dynamic_Viscosity(viscosity)
  call Control_Mod_Heat_Capacity(capacity)
  call Control_Mod_Mass_Density(density)
  call Control_Mod_Thermal_Conductivity(conductivity)

  end subroutine
