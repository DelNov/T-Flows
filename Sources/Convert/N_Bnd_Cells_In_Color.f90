!==============================================================================!
  integer function N_Bnd_Cells_In_Color(grid, bc, cell_data)
!------------------------------------------------------------------------------!
!   Counts and marks (with cell_data) boundary cells in given boundary color   !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: bc
  integer         :: cell_data(-grid % n_bnd_cells:grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, cnt
!==============================================================================!

  ! Nullify on entry
  cnt = 0
  cell_data(:) = 0

  do c = -grid % n_bnd_cells, -1
    if( grid % bnd_cond % color(c) .eq. bc ) then
      cnt = cnt + 1
      cell_data(c) = cell_data(c) + 1
    end if
  end do

  ! Return a sum of all marked boundary cells
  N_Bnd_Cells_In_Color = cnt

  end function

