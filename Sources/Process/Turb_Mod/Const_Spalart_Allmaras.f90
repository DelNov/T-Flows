!==============================================================================!
  subroutine Turb_Mod_Const_Spalart_Allmaras(turb)
!------------------------------------------------------------------------------!
!   Initializes constants for Spalar-Allmaras and DES turbulence models.       ! 
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Turb_Type),  target :: turb
!==============================================================================!

  turb % vis % sigma = TWO_THIRDS

  kappa  = 0.41
  c_b1   = 0.1355
  c_b2   = 0.622
  c_v1   = 7.1
  c_w1   = c_b1 / kappa**2 + (1 + c_b2) / turb % vis % sigma
  c_w2   = 0.3
  c_w3   = 2.0

  end subroutine
