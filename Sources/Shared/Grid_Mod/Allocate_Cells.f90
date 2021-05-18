!==============================================================================!
  subroutine Allocate_Cells(Grid, nc, nb)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type) :: Grid
  integer          :: nc    ! number of cells inside
  integer          :: nb    ! number of cells on the bounday
!-----------------------------------[Locals]-----------------------------------!
  integer :: c
!==============================================================================!

  ! Store number of cells and boundary cells
  Grid % n_cells     = nc
  Grid % n_bnd_cells = nb

  ! Allocate cell center coordinates and initialize to zero
  allocate(Grid % xc(-nb:nc));  Grid % xc(:) = 0.0
  allocate(Grid % yc(-nb:nc));  Grid % yc(:) = 0.0
  allocate(Grid % zc(-nb:nc));  Grid % zc(:) = 0.0

  ! Memory for cells' volumes, delta and wall distance
  allocate(Grid % vol           (-nb:nc));  Grid % vol           (:) = 0.0
  allocate(Grid % wall_dist     (-nb:nc));  Grid % wall_dist     (:) = 0.0
  allocate(Grid % cell_near_wall(-nb:nc));  Grid % cell_near_wall(:) = .false.

  ! Cells' nodes and neigboring cells
  allocate(Grid % cells_n(MAX_CELLS_N_NODES, -nb:nc));  Grid % cells_n(:,:) = 0
  allocate(Grid % cells_f(MAX_CELLS_N_FACES, -nb:nc));  Grid % cells_f(:,:) = 0
  allocate(Grid % cells_c(MAX_CELLS_N_CELLS, -nb:nc));  Grid % cells_c(:,:) = 0

  allocate(Grid % cells_bnd_face(-nb:-1));  Grid % cells_bnd_face(:) = 0

  ! Number of nodes, faces and cells at each cell
  ! (Actually, cells_n_faces and cells_n_cells should be the same)
  allocate(Grid % cells_n_nodes(-nb:nc));  Grid % cells_n_nodes(:) = 0
  allocate(Grid % cells_n_faces(-nb:nc));  Grid % cells_n_faces(:) = 0
  allocate(Grid % cells_n_cells(-nb:nc));  Grid % cells_n_cells(:) = 0

  ! Boundary condition color in a given direction
  ! (These go up to 6 because they are needed for
  !  non-polyhedral meshes creted in Gambit/Gmsh)
  allocate(Grid % cells_bnd_color(6, -nb:nc))

  ! Allocate processor i.d.
  allocate(Grid % comm % cell_proc(-nb:nc));  Grid % comm % cell_proc(:) = 0
  allocate(Grid % comm % cell_glo (-nb:nc))
  do c = -nb, nc
    Grid % comm % cell_glo(c) = c
  end do

  ! Allocate new and old numbers (this is so often used, maybe is better here)
  allocate(Grid % new_c(-nb:nc));  Grid % new_c(:) = 0
  allocate(Grid % old_c(-nb:nc));  Grid % old_c(:) = 0

  end subroutine
