!==============================================================================!
  subroutine Bnd_Cond_Ranges(Grid)
!------------------------------------------------------------------------------!
!   Allocates memory and finds the range (first and last boundary cell)        !
!   for each of the boundary condition colors.                                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type) :: Grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, color
!==============================================================================!

  ! Allocate memory
  allocate(Grid % bnd_cond % color_s_cell(Grid % n_bnd_cond))
  allocate(Grid % bnd_cond % color_e_cell(Grid % n_bnd_cond))

  ! Set non-realizable ranges
  Grid % bnd_cond % color_s_cell(:) = -1
  Grid % bnd_cond % color_e_cell(:) =  0

  ! Browse forward and backward to find first and last cell for each range
  do c = -1, -Grid % n_bnd_cells, -1
    color = Grid % bnd_cond % color(c)
    if(c < Grid % bnd_cond % color_e_cell(color)) then
      Grid % bnd_cond % color_e_cell(color) = c
    end if
  end do
  do c = -Grid % n_bnd_cells, -1
    color = Grid % bnd_cond % color(c)
    if(c > Grid % bnd_cond % color_e_cell(color)) then
      Grid % bnd_cond % color_s_cell(color) = c
    end if
  end do

  end subroutine
