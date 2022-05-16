!==============================================================================!
  subroutine Const_Reynolds_Stress(Turb)
!------------------------------------------------------------------------------!
!   Initializes constants for Reynolds Stress turbulence model.                ! 
!   (It is very simular to constants for Hanjalic-Jakirlic turbulence model    !
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
  c_l     =  0.161
  c_t     =  6.0
  c_nu    = 80.0
  g1      =  3.4
  g1_star =  1.8
  g2      =  4.2
  g3      =  0.8
  g3_star =  1.3
  g4      =  1.25
  g5      =  0.4

  Turb % kin % sigma = 1.0
  Turb % eps % sigma = 1.15
  Turb % uu  % sigma = 1.0
  Turb % vv  % sigma = 1.0
  Turb % ww  % sigma = 1.0
  Turb % uv  % sigma = 1.0
  Turb % uw  % sigma = 1.0
  Turb % vw  % sigma = 1.0

  end subroutine
