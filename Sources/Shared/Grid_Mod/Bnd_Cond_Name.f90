!==============================================================================!
  character(len=80) function Grid_Mod_Bnd_Cond_Name(grid, bnd_cell)
!------------------------------------------------------------------------------!
!   Provides a shortcut to obtain boundary condition type.                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: bnd_cell
!==============================================================================!

  Grid_Mod_Bnd_Cond_Name =  &
       grid % bnd_cond % name(grid % bnd_cond % color(bnd_cell))

  end function

