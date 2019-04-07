!==============================================================================!
  subroutine Turb_Mod_Alias_Stresses(turb, uu, vv, ww, uv, uw, vw)
!------------------------------------------------------------------------------!
!   Create aliases for Reynolds stress components.                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Turb_Type), target  :: turb
  type(Var_Type),  pointer :: uu, vv, ww, uv, uw, vw
!==============================================================================!

  uu => turb % uu
  vv => turb % vv
  ww => turb % ww
  uv => turb % uv
  uw => turb % uw
  vw => turb % vw

  end subroutine
