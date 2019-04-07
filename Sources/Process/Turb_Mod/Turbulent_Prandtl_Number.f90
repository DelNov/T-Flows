!==============================================================================!
  real function Turbulent_Prandtl_Number(grid, c)
!------------------------------------------------------------------------------!
!   Computes turbulent Prandtl number according to Kays, W. M. and             !
!   Crawford, M. E., Convective Heat and Mass Transfer,                        !
!   3rd edn. McGraw-Hill, New York, 1993.                                      !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod, only: TINY
  use Grid_Mod,  only: Grid_Type
  use Field_Mod, only: viscosity, conductivity, capacity 
  use Turb_Mod,  only: vis_t
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: c
  real            :: pr, pr_t_inf, pe_t
!==============================================================================!

  pr = viscosity * capacity / conductivity
  pr_t_inf = 0.85
  pe_t     = max(pr * vis_t(c) / viscosity, TINY)

  Turbulent_Prandtl_Number =                                            &
    1.0 / (   1.0/(2.0*pr_t_inf)                                        &
           +  0.3*pe_t * sqrt(1.0/pr_t_inf)                             &
           - (0.3*pe_t)**2 * (1.0-exp(-1.0/(0.3*pe_t*sqrt(pr_t_inf))))  &
          )

  end function
