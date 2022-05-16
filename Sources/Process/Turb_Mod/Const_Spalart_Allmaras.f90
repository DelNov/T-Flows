!==============================================================================!
  subroutine Const_Spalart_Allmaras(Turb)
!------------------------------------------------------------------------------!
!   Initializes constants for Spalar-Allmaras and DES turbulence models.       ! 
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type), target :: Turb
!==============================================================================!

  Turb % vis % sigma = TWO_THIRDS

  kappa  = 0.41
  c_b1   = 0.1355
  c_b2   = 0.622
  c_v1   = 7.1
  c_w1   = c_b1 / kappa**2 + (1 + c_b2) / Turb % vis % sigma
  c_w2   = 0.3
  c_w3   = 2.0

  end subroutine
