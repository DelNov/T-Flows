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

  c_mu   = 0.09
  c_mu25 = sqrt(sqrt(c_mu))
  c_mu75 = c_mu25**3
  kappa  = 0.41
  e_log  = 8.342

  ! Constants for buoyancy wall function
  c_mu_theta  = 0.1225
  c_mu_theta5 = 0.35
  kappa_theta = 0.38         ! von Karman constant for temperature

  ! Constants for GGDH turbulent flux
  c_theta = 0.22

  ! Constants for AFM turbulent flux
  afm_psi = 0.1
  afm_eta = 0.1

  end subroutine
