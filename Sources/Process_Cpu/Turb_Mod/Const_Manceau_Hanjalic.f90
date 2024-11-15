!==============================================================================!
  subroutine Const_Manceau_Hanjalic(Turb)
!------------------------------------------------------------------------------!
!   Initializes constants for Manceau Hanjalic Reynolds Stress turb. model     ! 
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
  Turb % c_l     =  0.161
  Turb % c_t     =  6.0
  Turb % c_nu    = 80.0
  Turb % g1      =  3.4
  Turb % g1_star =  1.8
  Turb % g2      =  4.2
  Turb % g3      =  0.8
  Turb % g3_star =  1.3
  Turb % g4      =  1.25
  Turb % g5      =  0.4

  ! Constants for GGDH turbulent flux
  Turb % c_theta =  0.22

  ! Constants for AFM turbulent flux
  Turb % afm_psi = 0.1
  Turb % afm_eta = 0.1

  Turb % kin % sigma = 1.0
  Turb % eps % sigma = 1.15
  Turb % uu % sigma  = 1.0
  Turb % vv % sigma  = 1.0
  Turb % ww % sigma  = 1.0
  Turb % uv % sigma  = 1.0
  Turb % uw % sigma  = 1.0
  Turb % vw % sigma  = 1.0

  end subroutine
