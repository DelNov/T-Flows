!==============================================================================!
  subroutine Grid_Mod_Bnd_Cond_Ranges(grid)
!------------------------------------------------------------------------------!
!   Allocates memory and finds the range (first and last boundary cell)        !
!   for each of the boundary condition colors.                                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, color
!==============================================================================!

  ! Allocate memory
  allocate(grid % bnd_cond % color_f(grid % n_bnd_cond))
  allocate(grid % bnd_cond % color_l(grid % n_bnd_cond))

  ! Set non-realizable ranges
  grid % bnd_cond % color_f(:) = -1
  grid % bnd_cond % color_l(:) =  0

  ! Browse forward and backward to find first and last cell for each range
  do c = -1, -grid % n_bnd_cells, -1
    color = grid % bnd_cond % color(c)
    if(c < grid % bnd_cond % color_l(color)) then
      grid % bnd_cond % color_l(color) = c
    end if
  end do
  do c = -grid % n_bnd_cells, -1
    color = grid % bnd_cond % color(c)
    if(c > grid % bnd_cond % color_l(color)) then
      grid % bnd_cond % color_f(color) = c
    end if
  end do

  end subroutine
