!==============================================================================!
  real function Turb_Mod_Prandtl_Number(turb, c)
!------------------------------------------------------------------------------!
!   Computes turbulent Prandtl number according to Kays, W. M. and             !
!   Crawford, M. E., Convective Heat and Mass Transfer,                        !
!   3rd edn. McGraw-Hill, New York, 1993.                                      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Turb_Type), target :: turb
  integer                 :: c
  real                    :: pr, pr_t_inf, pe_t
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: Flow
!==============================================================================!

  ! Take alias to the Flow to access its properties
  Flow => turb % pnt_flow

  pr = Flow % Prandtl_Number(c)
  pr_t_inf = 0.85
  pe_t     = max(pr * turb % vis_t(c) / Flow % viscosity(c), TINY)

  Turb_Mod_Prandtl_Number =                                             &
    1.0 / (   1.0/(2.0*pr_t_inf)                                        &
           +  0.3*pe_t * sqrt(1.0/pr_t_inf)                             &
           - (0.3*pe_t)**2 * (1.0-exp(-1.0/(0.3*pe_t*sqrt(pr_t_inf))))  &
          )

  end function
