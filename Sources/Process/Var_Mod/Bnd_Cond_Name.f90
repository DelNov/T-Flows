!==============================================================================!
  character(len=80) function Var_Mod_Bnd_Cond_Name(phi, bnd_cell)
!------------------------------------------------------------------------------!
!   Provides a shortcut to obtain boundary condition type.                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Var_Type) :: phi
  integer        :: bnd_cell
!------------------------------------------------------------------------------!

  Var_Mod_Bnd_Cond_Name = Grid_Mod_Bnd_Cond_Name(phi % pnt_grid, bnd_cell)

  end function

