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

  kappa = 0.41
  c_b1  = 0.1355
  c_b2  = 0.622
  c_v1  = 7.1
  c_w1  = c_b1 / kappa**2 + (1 + c_b2) / Turb % vis % sigma
  c_w2  = 0.3
  c_w3  = 2.0

  ! Constants for GGDH turbulent flux
  c_theta =  0.22

  ! Constants for AFM turbulent flux
  afm_psi = 0.1
  afm_eta = 0.1

  end subroutine
