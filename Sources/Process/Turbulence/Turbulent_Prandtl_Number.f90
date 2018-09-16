!==============================================================================!
  real function Turbulent_Prandtl_Number(grid, c)
!------------------------------------------------------------------------------!
!   Computes turbulent Prandtl number.                                         !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod,      only: TINY
  use Grid_Mod,       only: Grid_Type
  use Flow_Mod,       only: viscosity 
  use Turbulence_Mod, only: vis_t
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: c 
!==============================================================================!

  Turbulent_Prandtl_Number =                                                  &
    1.0 / (   0.5882                                                          &
            + 0.228  * (vis_t(c) / (viscosity + TINY))                        &
            - 0.0441 * (vis_t(c) / (viscosity + TINY))**2                     &
                     * (1.0 - exp(-5.165*( viscosity / (vis_t(c) + TINY) )))  &
          )

  end function
