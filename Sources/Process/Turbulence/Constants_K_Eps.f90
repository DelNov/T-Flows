!==============================================================================!
  subroutine Constants_K_Eps()
!------------------------------------------------------------------------------!
!   Initializes constants for k-eps turbulence model.                          ! 
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Rans_Mod
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  ! call Control_Mod_Turbulence_Model(.true.)
  ! call Control_Mod_Turbulence_Model_Variant(.true.)

  c_1e   = 1.44
  c_2e   = 1.92
  c_mu   = 0.09
  c_mu25 = sqrt(sqrt(c_mu))
  c_mu75 = c_mu25**3
  kappa  = 0.41
  e_log  = 8.342

  kin % sigma = 1.0
  eps % sigma = 1.3
! if(MODE .eq. LRe) then
!      c_1e = 1.55
!      c_2e = 2.0
! end if

  end subroutine
