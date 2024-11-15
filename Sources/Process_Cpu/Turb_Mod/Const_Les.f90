!==============================================================================!
  subroutine Const_Les(Turb)
!------------------------------------------------------------------------------!
!   Initializes constants for LES sub-grid scale models                        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type), target :: Turb
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Turb)
!==============================================================================!

  Turb % c_mu   = 0.09
  Turb % c_mu25 = sqrt(sqrt(Turb % c_mu))
  Turb % c_mu75 = Turb % c_mu25**3
  Turb % kappa  = 0.41
  Turb % e_log  = 8.342

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
