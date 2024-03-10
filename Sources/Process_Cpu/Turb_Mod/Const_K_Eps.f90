!==============================================================================!
  subroutine Const_K_Eps(Turb)
!------------------------------------------------------------------------------!
!   Initializes constants for k-eps turbulence model.                          ! 
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type), target :: Turb
!==============================================================================!

  c_1e    = 1.5
  c_2e    = 1.9
  c_mu    = 0.09
  c_mu25  = sqrt(sqrt(c_mu))
  c_mu75  = c_mu25**3
  kappa   = 0.41
  e_log   = 8.342

  ! Transport properties
  Turb % kin % sigma = 1.4
  Turb % eps % sigma = 1.4
  Turb % t2  % sigma = 1.0  ! copied from Const_K_Eps_Zeta_F

  c_mu_theta5 = 0.35
  kappa_theta = 0.38         ! von Karman constant for temperature

  ! Constants for GGDH turbulent flux
  c_theta = 0.2

  ! Constants for AFM turbulent flux
  afm_psi = 0.1
  afm_eta = 0.1

  end subroutine
