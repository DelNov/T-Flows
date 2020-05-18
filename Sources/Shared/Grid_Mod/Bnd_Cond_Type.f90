!==============================================================================!
  integer function Grid_Mod_Bnd_Cond_Type(grid, bnd_cell)
!------------------------------------------------------------------------------!
!   Provides a shortcut to obtain boundary condition type.                     !
!   Warning: avoid calling this function for temperature and species (scalars) !
!            since they may have different b.c. for the same color.            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: bnd_cell
!==============================================================================!

  Grid_Mod_Bnd_Cond_Type =  &
       grid % bnd_cond % type(grid % bnd_cond % color(bnd_cell))

  end function

