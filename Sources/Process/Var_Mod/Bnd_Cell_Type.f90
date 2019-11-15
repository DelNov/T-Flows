!==============================================================================!
  integer function Var_Mod_Bnd_Cell_Type(phi, bnd_cell)
!------------------------------------------------------------------------------!
!   Provides a shortcut to obtain boundary condition type.                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Var_Type) :: phi
  integer        :: bnd_cell
!------------------------------------------------------------------------------!

  Var_Mod_Bnd_Cell_Type = phi % bnd_cell_type(bnd_cell)

  end function

