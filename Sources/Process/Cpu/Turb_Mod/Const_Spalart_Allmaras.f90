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

  Turb % kappa = 0.41
  Turb % e_log = 8.34
  Turb % c_b1  = 0.1355
  Turb % c_b2  = 0.622
  Turb % c_v1  = 7.1
  Turb % c_w1  = Turb % c_b1 / Turb % kappa**2  &
               + (1 + Turb % c_b2) / Turb % vis % sigma
  Turb % c_w2  = 0.3
  Turb % c_w3  = 2.0

  Turb % c_mu_theta5 = 0.35
  Turb % kappa_theta = 0.38  ! von Karman constant for temperature

  ! Constants for GGDH turbulent flux
  Turb % c_theta =  0.22

  ! Constants for AFM turbulent flux
  Turb % afm_psi = 0.1
  Turb % afm_eta = 0.1

  end subroutine
