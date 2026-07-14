!==============================================================================!
  integer function N_Bnd_Cells_In_Region(Convert, Grid, bc, bnd_cell_flag)
!------------------------------------------------------------------------------!
!>  Counts and marks boundary cells in a specified boundary condition (bc)
!>  category within a grid.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Convert_Type) :: Convert  !! parent class
  type(Grid_Type)     :: Grid     !! grid being converted
  integer             :: bc       !! boundary condition rank (number)
  integer             :: bnd_cell_flag(-Grid % n_bnd_cells  &
                                       :Grid % n_cells)  !! stored flag on
                                                         !! boundary cells
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, cnt
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Convert)
!==============================================================================!

  ! Nullify on entry
  bnd_cell_flag(:) = 0

  cnt = 0
  do c = -Grid % n_bnd_cells, -1
    if( Grid % region % at_cell(c) .eq. bc ) then
      cnt = cnt + 1
      bnd_cell_flag(c) = bnd_cell_flag(c) + 1
    end if
  end do

  ! Return a sum of all marked boundary cells
  N_Bnd_Cells_In_Region = cnt

  end function

