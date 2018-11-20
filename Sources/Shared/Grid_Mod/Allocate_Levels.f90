!==============================================================================!
  subroutine Grid_Mod_Allocate_Levels(grid)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer         :: lev
!==============================================================================!

  ! At each level, allocate memory for individual data
  do lev = 0, MAX_MG_LEV
    allocate(grid % level(lev) % cell(grid % n_cells))
    allocate(grid % level(lev) % face(grid % n_faces))
    allocate(grid % level(lev) % faces_c(2, grid % n_faces))
    allocate(grid % level(lev) % n_finest_cells(grid % n_cells))
    allocate(grid % level(lev) % xc(grid % n_cells))
    allocate(grid % level(lev) % yc(grid % n_cells))
    allocate(grid % level(lev) % zc(grid % n_cells))
  end do

  end subroutine
