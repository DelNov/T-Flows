!==============================================================================!
  subroutine Turb_Mod_Const_Hanjalic_Jakirlic(turb)
!------------------------------------------------------------------------------!
!   Initializes constants for Hanjalic-Jakirlic turbulence model.              ! 
!   (It is very simular to constants for Reynolds Stress turbulence model.     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Turb_Type), target :: turb
!==============================================================================!

  c_1e    =  1.44
  c_2e    =  1.83
  c_3e    =  0.55
  c_mu    =  0.09
  c_mu_d  =  0.21
  c_mu25  = sqrt(sqrt(c_mu))
  c_mu75  = c_mu25**3
  c_theta = 0.25

  turb % kin % sigma = 1.0
  turb % eps % sigma = 1.0
  turb % uu  % sigma = 1.0
  turb % vv  % sigma = 1.0
  turb % ww  % sigma = 1.0
  turb % uv  % sigma = 1.0
  turb % uw  % sigma = 1.0
  turb % vw  % sigma = 1.0

  end subroutine
