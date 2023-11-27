!==============================================================================!
  character(SL) function Var_Mod_Bnd_Cond_Name(phi, bnd_cell)
!------------------------------------------------------------------------------!
!>  This function provides a convenient and quick way to retrieve the name
!>  (in the form of a string) of the boundary condition associated with a
!>  specific boundary cell.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Var_Type) :: phi       !! variable object under consideration
  integer        :: bnd_cell  !! boundary cell number
!------------------------------------------------------------------------------!

  Var_Mod_Bnd_Cond_Name = phi % pnt_grid % Bnd_Cond_Name_At_Cell(bnd_cell)

  end function

