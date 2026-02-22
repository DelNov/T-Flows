!==============================================================================!
  subroutine Const_K_Omega_Sst(Turb)
!------------------------------------------------------------------------------!
!   Initializes constants for k-omega SST turbulence model.
!
!   NOTE:
!   - This routine intentionally keeps the existing wall-function and
!     heat-transfer/scalar-flux coefficients used elsewhere in the code
!     (kappa, e_log, c_mu_theta5, kappa_theta, c_theta, afm_psi, afm_eta, ...).
!   - SST-specific constants added here are the minimal set needed by the
!     current k-omega-SST implementation (k/omega sources + vis_t limiter).
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type), target :: Turb
!==============================================================================!

  !---------------------------!
  !   SST model constants     !
  !---------------------------!
  ! Menter SST (common baseline values)
  Turb % beta_star = 0.09
  Turb % a1        = 0.31

  ! For the current implementation we keep scalar beta (used in omega sink term)
  ! and beta_1 (used for wall omega in integrate-to-wall / wall-function blend).
  Turb % beta1    = 0.075
  Turb % beta2    = 0.0828

  Turb % gamma1 = 5.0 / 9.0     ! ≈ 0.5555556
  Turb % gamma2 = 0.44

  ! If you later implement full SST blending (F1/F2), these become cell-wise:
  !   beta  = F1*beta1 + (1-F1)*beta2
  !   gamma = F1*gamma1 + (1-F1)*gamma2
  !   sigma_k, sigma_w similarly.
  !
  ! For now, diffusion is handled through phi%sigma in Compute_Variable, so we
  ! keep the same (global) sigmas as in your original file.

  !-------------------------------!
  !   Wall-function coefficients  !
  !-------------------------------!
  Turb % c_mu   = 0.09
  Turb % c_mu25 = sqrt(sqrt(Turb % c_mu))
  Turb % c_mu75 = Turb % c_mu25**3
  Turb % kappa  = 0.41
  Turb % e_log  = 8.342

  !------------------------!
  !   Transport properties !
  !------------------------!
  Turb % sig_k1 = 0.85
  Turb % sig_k2 = 1.00
  Turb % sig_w1 = 0.50
  Turb % sig_w2 = 0.856
  Turb % t2 % sigma = 1.0  ! copied from Const_K_Eps_Zeta_F

  !----------------------------------------!
  !   Heat transfer / scalar-flux models   !
  !----------------------------------------!
  Turb % c_mu_theta5 = 0.35
  Turb % kappa_theta = 0.38   ! von Karman constant for temperature

  ! Constants for GGDH turbulent flux
  Turb % c_theta = 0.2

  ! Constants for AFM turbulent flux
  Turb % afm_psi = 0.1
  Turb % afm_eta = 0.1

  end subroutine Const_K_Omega_Sst
