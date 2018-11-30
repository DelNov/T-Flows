!==============================================================================!
  subroutine Grid_Mod_Allocate_Levels(grid)
!------------------------------------------------------------------------------!
!   This function is called from Process                                       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: lev, n_cells, n_faces, c, s
!==============================================================================!

  !--------------------------------!
  !                                !
  !   Allocate memory for levels   !
  !                                !
  !--------------------------------!

  ! At each level, allocate memory for individual data
  do lev = 1, grid % n_levels

    ! Tentative number of cells and faces
    n_cells = grid % level(lev) % n_cells
    n_faces = grid % level(lev) % n_faces

    ! Connectivity with finest grid
    allocate(grid % level(lev) % cell(grid % n_cells))
    allocate(grid % level(lev) % face(grid % n_faces))
    grid % level(lev) % cell(:) = 0
    grid % level(lev) % face(:) = 0

    ! Connectivity with coarser grid
    allocate(grid % level(lev) % coarser_c(n_cells))
    grid % level(lev) % coarser_c(:) = 0

    ! Intra-level connectivity and geometry
    allocate(grid % level(lev) % faces_c(2, n_faces))
    allocate(grid % level(lev) % n_finest_cells(n_cells))
    allocate(grid % level(lev) % xc     (   n_cells))
    allocate(grid % level(lev) % yc     (   n_cells))
    allocate(grid % level(lev) % zc     (   n_cells))
  end do

  end subroutine
