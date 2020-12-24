!==============================================================================!
  subroutine Grid_Mod_Allocate_Cells(grid, nc, nb)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: nc        ! number of cells inside
  integer         :: nb        ! number of cells on the bounday
!-----------------------------------[Locals]-----------------------------------!
  integer :: c
!==============================================================================!

  ! Store number of cells and boundary cells
  grid % n_cells     = nc
  grid % n_bnd_cells = nb

  ! Allocate cell center coordinates and initialize to zero
  allocate(grid % xc(-nb:nc));  grid % xc(:) = 0.0
  allocate(grid % yc(-nb:nc));  grid % yc(:) = 0.0
  allocate(grid % zc(-nb:nc));  grid % zc(:) = 0.0

  ! Memory for cells' volumes, delta and wall distance
  allocate(grid % vol           (-nb:nc));  grid % vol           (:) = 0.0
  allocate(grid % wall_dist     (-nb:nc));  grid % wall_dist     (:) = 0.0
  allocate(grid % cell_near_wall(-nb:nc));  grid % cell_near_wall(:) = .false.

  ! Cells' nodes and neigboring cells
  allocate(grid % cells_n(MAX_CELLS_N_NODES, -nb:nc));  grid % cells_n(:,:) = 0
  allocate(grid % cells_f(MAX_CELLS_N_FACES, -nb:nc));  grid % cells_f(:,:) = 0
  allocate(grid % cells_c(MAX_CELLS_N_CELLS, -nb:nc));  grid % cells_c(:,:) = 0

  allocate(grid % cells_bnd_face(-nb:-1));  grid % cells_bnd_face(:) = 0

  ! Number of nodes, faces and cells at each cell
  ! (Actually, cells_n_faces and cells_n_cells should be the same)
  allocate(grid % cells_n_nodes(-nb:nc));  grid % cells_n_nodes(:) = 0
  allocate(grid % cells_n_faces(-nb:nc));  grid % cells_n_faces(:) = 0
  allocate(grid % cells_n_cells(-nb:nc));  grid % cells_n_cells(:) = 0

  ! Boundary condition color in a given direction
  ! (These go up to 6 because they are needed for
  !  non-polyhedral meshes creted in Gambit/Gmsh)
  allocate(grid % cells_bnd_color(6, -nb:nc))

  ! Allocate processor i.d.
  allocate(grid % comm % cell_proc(-nb:nc));  grid % comm % cell_proc(:) = 0
  allocate(grid % comm % cell_glo (-nb:nc))
  do c = -nb, nc
    grid % comm % cell_glo(c) = c
  end do

  ! Allocate new and old numbers (this is so often used, maybe is better here)
  allocate(grid % new_c(-nb:nc));  grid % new_c(:) = 0
  allocate(grid % old_c(-nb:nc));  grid % old_c(:) = 0

  end subroutine
