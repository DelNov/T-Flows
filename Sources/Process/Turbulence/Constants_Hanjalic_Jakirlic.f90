!==============================================================================!
  subroutine Constants_Hanjalic_Jakirlic()
!------------------------------------------------------------------------------!
!   Initializes constants for Hanjalic-Jakirlic turbulence model.              ! 
!   (It is very simular to constants for Reynolds Stress turbulence model.     !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Rans_Mod
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  ! call Control_Mod_Turbulence_Model(.true.)
  ! call Control_Mod_Turbulence_Model_Variant(.true.)

  c_1e    =  1.44
  c_2e    =  1.83
  c_3e    =  0.55
  c_mu    =  0.09
  c_mu_d  =  0.21
  c_mu25  = sqrt(sqrt(c_mu))
  c_mu75  = c_mu25**3

  kin % sigma = 1.0
  eps % sigma = 1.0 
  uu  % sigma = 1.0
  vv  % sigma = 1.0
  ww  % sigma = 1.0
  uv  % sigma = 1.0
  uw  % sigma = 1.0
  vw  % sigma = 1.0

  c_theta = 0.25

  end subroutine
