!==============================================================================!
  integer function Grid_Mod_Bnd_Cond_Type(grid, bnd_cell)
!------------------------------------------------------------------------------!
!   Provides a shortcut to obtain boundary condition type.                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: bnd_cell
!==============================================================================!

  Grid_Mod_Bnd_Cond_Type =  &
       grid % bnd_cond % type(grid % bnd_cond % color(bnd_cell))

  end function

