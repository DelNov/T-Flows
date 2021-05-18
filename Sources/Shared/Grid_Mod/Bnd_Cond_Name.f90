!==============================================================================!
  character(SL) function Bnd_Cond_Name(Grid, bnd_cell)
!------------------------------------------------------------------------------!
!   Provides a shortcut to obtain boundary condition type.                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type) :: Grid
  integer          :: bnd_cell
!==============================================================================!

  Bnd_Cond_Name =  &
       Grid % bnd_cond % name(Grid % bnd_cond % color(bnd_cell))

  end function

