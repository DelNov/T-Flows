!==============================================================================!
  subroutine Const_Hanjalic_Jakirlic(Turb)
!------------------------------------------------------------------------------!
!   Initializes constants for Hanjalic-Jakirlic turbulence model.              ! 
!   (It is very simular to constants for Reynolds Stress turbulence model.     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type), target :: Turb
!==============================================================================!

  Turb % c_1e    =  1.44
  Turb % c_2e    =  1.83
  Turb % c_3e    =  0.55
  Turb % c_mu    =  0.09
  Turb % c_mu_d  =  0.21
  Turb % c_mu25  = sqrt(sqrt(Turb % c_mu))
  Turb % c_mu75  = Turb % c_mu25**3

  Turb % kin % sigma = 1.0
  Turb % eps % sigma = 1.0
  Turb % uu  % sigma = 1.0
  Turb % vv  % sigma = 1.0
  Turb % ww  % sigma = 1.0
  Turb % uv  % sigma = 1.0
  Turb % uw  % sigma = 1.0
  Turb % vw  % sigma = 1.0

  ! Constants for GGDH turbulent flux
  Turb % c_theta =  0.255555

  ! Constants for AFM turbulent flux
  Turb % afm_psi = 0.1
  Turb % afm_eta = 0.1

  end subroutine
