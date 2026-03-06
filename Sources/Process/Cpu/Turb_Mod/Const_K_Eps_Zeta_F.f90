!==============================================================================!
  subroutine Const_K_Eps_Zeta_F(Turb)
!------------------------------------------------------------------------------!
!   Initializes constants for k-eps-zeta-f and hybrid k-eps-zeta-f             !
!   turbulence models.                                                         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type), target :: Turb
!==============================================================================!

  Turb % c_1e    =  1.4
  Turb % c_2e    =  1.9
  Turb % c_mu    =  0.09
  Turb % c_mu_d  =  0.22
  Turb % c_mu25  = sqrt(sqrt(Turb % c_mu))
  Turb % c_mu75  = Turb % c_mu25**3
  Turb % kappa   =  0.41
  Turb % e_log   =  8.342
  Turb % c_l     =  0.36
  Turb % c_t     =  6.0
  Turb % c_nu    = 85.0
  Turb % alpha   =  0.012
  Turb % c_f1    =  1.4
  Turb % c_f2    =  0.3

  ! Transport properties
  Turb % kin  % sigma = 1.0
  Turb % eps  % sigma = 1.3
  Turb % zeta % sigma = 1.2
  Turb % t2   % sigma = 1.0

  ! Constants for buoyancy wall function
  Turb % c_mu_theta  = 0.1225
  Turb % c_mu_theta5 = 0.35
  Turb % kappa_theta = 0.38         ! von Karman constant for temperature

  ! Constants for GGDH turbulent flux
  Turb % c_theta = 0.22

  ! Constants for AFM turbulent flux
  Turb % afm_psi = 0.1
  Turb % afm_eta = 0.1

  end subroutine
