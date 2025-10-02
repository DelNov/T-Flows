!==============================================================================!
  subroutine Alias_K_Eps_Zeta_F(Turb, kin, eps, zeta, f22)
!------------------------------------------------------------------------------!
!   Create aliases for k, epsilon, zeta and f22                                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type), target  :: Turb
  type(Var_Type),   pointer :: kin, eps, zeta, f22
!==============================================================================!

  kin  => Turb % kin
  eps  => Turb % eps
  zeta => Turb % zeta
  f22  => Turb % f22

  end subroutine
