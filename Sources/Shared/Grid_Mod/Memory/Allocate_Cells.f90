!==============================================================================!
  subroutine Allocate_Cells(Grid, nc, nb)
!------------------------------------------------------------------------------!
!>  Allocates memory for cell-based data (arrays and matrices), for geometrical
!>  (xc, yc, zc, vol ...) and connectivity data (cells_n, cells_f, ...).
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type)    :: Grid  !! grid under consideration
  integer, intent(in) :: nc    !! number of cells inside
  integer, intent(in) :: nb    !! number of cells on the bounday
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

  ! Memory for cells' inertia tensors
  allocate(Grid % ixx(-nb:nc));  Grid % ixx(:) = 0.0
  allocate(Grid % iyy(-nb:nc));  Grid % iyy(:) = 0.0
  allocate(Grid % izz(-nb:nc));  Grid % izz(:) = 0.0
  allocate(Grid % ixy(-nb:nc));  Grid % ixy(:) = 0.0
  allocate(Grid % ixz(-nb:nc));  Grid % ixz(:) = 0.0
  allocate(Grid % iyz(-nb:nc));  Grid % iyz(:) = 0.0

  ! Allocate as litle as possible
  allocate(Grid % cells_n(4, -nb:nc));  Grid % cells_n(:,:) = 0
  allocate(Grid % cells_f(4, -nb:nc));  Grid % cells_f(:,:) = 0
  if(PROGRAM_NAME .eq. "Generate" .or.  &
     PROGRAM_NAME .eq. "Convert") then
    allocate(Grid % cells_c(4, -nb:nc));  Grid % cells_c(:,:) = 0
  end if

  allocate(Grid % cells_bnd_face(-nb:-1));  Grid % cells_bnd_face(:) = 0

  ! Number of nodes, faces and cells at each cell
  ! (Actually, cells_n_faces and cells_n_cells should be the same)
  allocate(Grid % cells_n_nodes(-nb:nc));  Grid % cells_n_nodes(:) = 0
  allocate(Grid % cells_n_faces(-nb:nc));  Grid % cells_n_faces(:) = 0
  if(PROGRAM_NAME .eq. "Generate" .or.  &
     PROGRAM_NAME .eq. "Convert") then
    allocate(Grid % cells_n_cells(-nb:nc));  Grid % cells_n_cells(:) = 0
  end if

  ! Boundary condition region in a given direction
  ! (These go up to 6 because they are needed for
  !  non-polyhedral meshes creted in Gambit/Gmsh)
  allocate(Grid % cells_bnd_region(6, -nb:nc))

  ! Allocate cell-based porous regions
  allocate(Grid % por(-nb:nc));  Grid % por(:) = 0

  ! Allocate processor i.d.
  allocate(Grid % Comm % cell_proc(-nb:nc));  Grid % Comm % cell_proc(:) = 0
  allocate(Grid % Comm % cell_glo (-nb:nc))
  do c = -nb, nc
    Grid % Comm % cell_glo(c) = c
  end do

  ! Allocate thread i.d.
  allocate(Grid % Omp % cell_thread(-nb:nc));  Grid % Omp % cell_thread(:) = 0

  ! Allocate new and old numbers (this is so often used, maybe is better here)
  if(PROGRAM_NAME .ne. "Process") then
    allocate(Grid % new_c(-nb:nc));  Grid % new_c(:) = 0
    allocate(Grid % old_c(-nb:nc));  Grid % old_c(:) = 0
  end if

  end subroutine
