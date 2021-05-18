!==============================================================================!
  integer function Bnd_Cond_Type(Grid, bnd_cell)
!------------------------------------------------------------------------------!
!   Provides a shortcut to obtain boundary condition type.                     !
!   Warning: avoid calling this function for temperature and species (scalars) !
!            since they may have different b.c. for the same color.            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type) :: Grid
  integer          :: bnd_cell
!==============================================================================!

  Bnd_Cond_Type =  &
       Grid % bnd_cond % type(Grid % bnd_cond % color(bnd_cell))

  end function

