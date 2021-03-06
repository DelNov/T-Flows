!==============================================================================!
  subroutine Turb_Mod_Alias_Heat_Fluxes(turb, ut, vt, wt)
!------------------------------------------------------------------------------!
!   Create aliases for turbulent heat fluxes.                                  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Turb_Type), target  :: turb
  type(Var_Type),  pointer :: ut, vt, wt
!==============================================================================!

  ut => turb % ut
  vt => turb % vt
  wt => turb % wt

  end subroutine
