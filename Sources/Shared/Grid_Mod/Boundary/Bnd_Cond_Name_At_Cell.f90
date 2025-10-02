!==============================================================================!
  character(SL) function Bnd_Cond_Name_At_Cell(Grid, cell)
!------------------------------------------------------------------------------!
!>  Provides a shortcut to obtain boundary condition name at boundary cell.
!>  It is a bit of the relict of the past, it would be better to browse
!>  through bundary regions first, check the name of the region, and then
!>  browse through all the cells in that region.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type)    :: Grid  !! grid under consideration
  integer, intent(in) :: cell  !! cell number
!==============================================================================!

  Bnd_Cond_Name_At_Cell = Grid % region % name(Grid % region % at_cell(cell))

  end function

