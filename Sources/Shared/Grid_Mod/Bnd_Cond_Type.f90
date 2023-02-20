!==============================================================================!
  integer function Bnd_Cond_Type(Grid, cell)
!------------------------------------------------------------------------------!
!   Provides a shortcut to obtain boundary condition type.                     !
!   Warning: avoid calling this function for temperature and species (scalars) !
!            since they may have different b.c. for the same region.           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type)    :: Grid
  integer, intent(in) :: cell
!==============================================================================!

  Bnd_Cond_Type = Grid % region % type(Grid % region % at_cell(cell))

  end function

