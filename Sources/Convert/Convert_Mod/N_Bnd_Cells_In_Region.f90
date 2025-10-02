!==============================================================================!
  integer function N_Bnd_Cells_In_Region(Convert, Grid, bc, cell_data)
!------------------------------------------------------------------------------!
!>  Counts and marks boundary cells in a specified boundary condition (bc)
!>  category within a grid.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Convert_Type) :: Convert  !! parent class
  type(Grid_Type)     :: Grid     !! grid being converted
  integer             :: bc       !! boundary condition rank (number)
  integer             :: cell_data(-Grid % n_bnd_cells  &
                                   :Grid % n_cells)  !! stored data on
                                                     !! boundary cells
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, cnt
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Convert)
!==============================================================================!

  ! Nullify on entry
  cnt = 0
  cell_data(:) = 0

  do c = -Grid % n_bnd_cells, -1
    if( Grid % region % at_cell(c) .eq. bc ) then
      cnt = cnt + 1
      cell_data(c) = cell_data(c) + 1
    end if
  end do

  ! Return a sum of all marked boundary cells
  N_Bnd_Cells_In_Region = cnt

  end function

