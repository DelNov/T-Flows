!==============================================================================!
  character(SL) function Bnd_Cond_Name_At_Cell(Grid, cell)
!------------------------------------------------------------------------------!
!   Provides a shortcut to obtain boundary condition name at boundary cell.    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type)    :: Grid
  integer, intent(in) :: cell
!==============================================================================!

  Bnd_Cond_Name_At_Cell = Grid % region % name(Grid % region % at_cell(cell))

  end function

