!==============================================================================!
  character(SL) function Bnd_Cond_Name(Grid, cell)
!------------------------------------------------------------------------------!
!   Provides a shortcut to obtain boundary condition type.                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type)    :: Grid
  integer, intent(in) :: cell
!==============================================================================!

  Bnd_Cond_Name = Grid % region % name(Grid % region % at_cell(cell))

  end function

