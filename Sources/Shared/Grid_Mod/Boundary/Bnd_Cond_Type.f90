!==============================================================================!
  integer function Bnd_Cond_Type(Grid, cell)
!------------------------------------------------------------------------------!
!>  Provides a shortcut to obtain boundary condition type for a cell.
!>  It is a bit of the relict of the past, it would be better to browse
!>  through bundary regions first, check the type of the region, and then
!>  browse through all the cells in that region.
!   Warning: avoid calling this function for temperature and species (scalars) !
!            since they may have different b.c. for the same region.           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type)    :: Grid  !! grid under consideration
  integer, intent(in) :: cell  !! cell number
!==============================================================================!

  Bnd_Cond_Type = Grid % region % type(Grid % region % at_cell(cell))

  end function

