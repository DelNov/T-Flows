!==============================================================================!
  subroutine Multiphase_Mod_Alias_Vof(mult, vof)
!------------------------------------------------------------------------------!
!   Creates aliases for volume fraction.  (Actually, only one alias.)          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target  :: mult
  type(Var_Type),        pointer :: vof
!==============================================================================!

  vof => mult % vof

  end subroutine
