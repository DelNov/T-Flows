!==============================================================================!
  subroutine Turb_Mod_Alias_K_Eps_Zeta_F(turb, kin, eps, zeta, f22)
!------------------------------------------------------------------------------!
!   Create aliases for k, epsilon, zeta and f22                                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Turb_Type), target  :: turb
  type(Var_Type),  pointer :: kin, eps, zeta, f22
!==============================================================================!

  kin  => turb % kin
  eps  => turb % eps
  zeta => turb % zeta
  f22  => turb % f22

  end subroutine
