!==============================================================================!
  integer function Var_Mod_Bnd_Cond_Type(phi, bnd_cell)
!------------------------------------------------------------------------------!
!>  Similar to Var_Mod_Bnd_Cond_Name, this function provides a direct way to
!>  obtain the type of boundary condition assigned to a particular boundary
!>  cell. It focuses on the type (numerical representation) rather than a
!>  These numerical representations of types are defined in Region_Mod.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Var_Type) :: phi       !! variable object under consideration
  integer        :: bnd_cell  !! boundary cell number
!------------------------------------------------------------------------------!

  Var_Mod_Bnd_Cond_Type = phi % bnd_cond_type(bnd_cell)

  end function

