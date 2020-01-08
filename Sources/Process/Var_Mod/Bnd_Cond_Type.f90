!==============================================================================!
  integer function Var_Mod_Bnd_Cond_Type(phi, bnd_cell)
!------------------------------------------------------------------------------!
!   Provides a shortcut to obtain boundary condition type.                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Var_Type) :: phi
  integer        :: bnd_cell
!------------------------------------------------------------------------------!

  Var_Mod_Bnd_Cond_Type = phi % bnd_cond_type(bnd_cell)

  end function

