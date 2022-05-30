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

  c_1e    =  1.44
  c_2e    =  1.83
  c_3e    =  0.55
  c_mu    =  0.09
  c_mu_d  =  0.21
  c_mu25  = sqrt(sqrt(c_mu))
  c_mu75  = c_mu25**3

  Turb % kin % sigma = 1.0
  Turb % eps % sigma = 1.0
  Turb % uu  % sigma = 1.0
  Turb % vv  % sigma = 1.0
  Turb % ww  % sigma = 1.0
  Turb % uv  % sigma = 1.0
  Turb % uw  % sigma = 1.0
  Turb % vw  % sigma = 1.0

  ! Constants for GGDH turbulent flux
  c_theta =  0.255555

  ! Constants for AFM turbulent flux
  afm_psi      = 0.1 
  afm_eta      = 0.1 

  end subroutine
