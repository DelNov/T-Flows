!==============================================================================!
  character(SL) function Var_Mod_Bnd_Cond_Name(phi, bnd_cell)
!------------------------------------------------------------------------------!
!   Provides a shortcut to obtain boundary condition type.                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Var_Type) :: phi
  integer        :: bnd_cell
!------------------------------------------------------------------------------!

  Var_Mod_Bnd_Cond_Name = phi % pnt_grid % Bnd_Cond_Name_At_Cell(bnd_cell)

  end function

