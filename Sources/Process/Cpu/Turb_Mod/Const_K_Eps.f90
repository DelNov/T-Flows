!==============================================================================!
  subroutine Const_K_Eps(Turb)
!------------------------------------------------------------------------------!
!   Initializes constants for k-eps turbulence model.                          ! 
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type), target :: Turb
!==============================================================================!

  Turb % c_1e    = 1.5
  Turb % c_2e    = 1.9
  Turb % c_mu    = 0.09
  Turb % c_mu25  = sqrt(sqrt(Turb % c_mu))
  Turb % c_mu75  = Turb % c_mu25**3
  Turb % kappa   = 0.41
  Turb % e_log   = 8.342

  ! Transport properties
  Turb % kin % sigma = 1.4
  Turb % eps % sigma = 1.4
  Turb % t2  % sigma = 1.0  ! copied from Const_K_Eps_Zeta_F

  Turb % c_mu_theta5 = 0.35
  Turb % kappa_theta = 0.38         ! von Karman constant for temperature

  ! Constants for GGDH turbulent flux
  Turb % c_theta = 0.2

  ! Constants for AFM turbulent flux
  Turb % afm_psi = 0.1
  Turb % afm_eta = 0.1

  end subroutine
