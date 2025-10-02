!==============================================================================!
  real function Prandtl_Turb(Turb, c)
!------------------------------------------------------------------------------!
!   Computes turbulent Prandtl number according to Kays, W. M. and             !
!   Crawford, M. E., Convective Heat and Mass Transfer, 3rd ed. McGraw-Hill,   !
!   New York, 1993 if it is not an LES, otherwise return simply 0.4.           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type), target :: Turb
  integer,      intent(in) :: c
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: Flow
  real                      :: pr, pr_t_inf, pe_t
!==============================================================================!

  ! Take alias to the Flow to access its properties
  Flow => Turb % pnt_flow

  ! In case an LES model is used, set it simply to 0.4
  if(Turb % Les()) then
    Prandtl_Turb = 0.4

  ! Otherwise, use the correlation proposed by Kays and Crawford
  else
    pr = Flow % Prandtl_Numb(c)
    pr_t_inf = 0.85
    pe_t     = max(pr * Turb % vis_t(c) / Flow % viscosity(c), TINY)

    Prandtl_Turb =                                                        &
      1.0 / (   1.0/(2.0*pr_t_inf)                                        &
             +  0.3*pe_t * sqrt(1.0/pr_t_inf)                             &
             - (0.3*pe_t)**2 * (1.0-exp(-1.0/(0.3*pe_t*sqrt(pr_t_inf))))  &
            )
  end if

  end function
