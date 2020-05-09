!==============================================================================!
  subroutine Grid_Mod_Allocate_Cells(grid, nc, nb)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: nc        ! number of cells inside
  integer         :: nb        ! number of cells on the bounday
!==============================================================================!

  ! Store number of cells and boundary cells
  grid % n_cells     = nc
  grid % n_bnd_cells = nb

  ! Allocate cell center coordinates and initialize to zero
  allocate(grid % xc(-nb:nc));  grid % xc = 0.0
  allocate(grid % yc(-nb:nc));  grid % yc = 0.0
  allocate(grid % zc(-nb:nc));  grid % zc = 0.0

  ! Memory for cells' volumes, delta and wall distance
  allocate(grid % vol           (-nb:nc));  grid % vol            = 0.0
  allocate(grid % wall_dist     (-nb:nc));  grid % wall_dist      = 0.0
  allocate(grid % cell_near_wall(-nb:nc));  grid % cell_near_wall = .false.

  ! Cells' nodes and neigboring cells
  allocate(grid % cells_n( 8, -nb:nc));  grid % cells_n = 0
  allocate(grid % cells_f( 6, -nb:nc));  grid % cells_f = 0
  allocate(grid % cells_c(24, -nb:nc));  grid % cells_c = 0

  allocate(grid % cells_bnd_face(-nb:-1));  grid % cells_bnd_face = 0

  ! Number of nodes at each cell (determines cell's shape really)
  allocate(grid % cells_n_nodes(-nb:nc));  grid % cells_n_nodes = 0

  ! Number of faces at each cell
  allocate(grid % cells_n_faces(-nb:nc));  grid % cells_n_faces = 0

  ! Boundary condition color in a given direction
  allocate(grid % cells_bnd_color(6, -nb:nc))

  ! Allocate processor i.d.
  allocate(grid % comm % cell_proc(-nb:nc));  grid % comm % cell_proc(:) = 0

  end subroutine
