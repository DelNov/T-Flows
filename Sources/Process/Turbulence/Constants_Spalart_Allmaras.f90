!==============================================================================!
  subroutine Constants_Spalart_Allmaras()
!------------------------------------------------------------------------------!
!   Initializes constants for Spalar-Allmaras and DES turbulence models.       ! 
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Rans_Mod
  use Control_Mod
  use Const_Mod, only: TWO_THIRDS
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  ! call Control_Mod_Turbulence_Model(.true.)
  ! call Control_Mod_Turbulence_Model_Variant(.true.)

  vis % sigma = TWO_THIRDS

  kappa  = 0.41
  c_b1   = 0.1355
  c_b2   = 0.622
  c_v1   = 7.1
  c_w1   = c_b1 / kappa**2 + (1 + c_b2) / vis % sigma
  c_w2   = 0.3
  c_w3   = 2.0

  end subroutine
