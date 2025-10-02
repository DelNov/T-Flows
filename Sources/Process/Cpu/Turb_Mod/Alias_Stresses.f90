!==============================================================================!
  subroutine Alias_Stresses(Turb, uu, vv, ww, uv, uw, vw)
!------------------------------------------------------------------------------!
!   Create aliases for Reynolds stress components.                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type), target  :: Turb
  type(Var_Type),   pointer :: uu, vv, ww, uv, uw, vw
!==============================================================================!

  uu => Turb % uu
  vv => Turb % vv
  ww => Turb % ww
  uv => Turb % uv
  uw => Turb % uw
  vw => Turb % vw

  end subroutine
