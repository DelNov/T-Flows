!==============================================================================!
  subroutine Grid_Mod_Allocate_Faces(grid, nf, ns)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: nf    ! number of faces in the grid   
  integer         :: ns    ! number of shadow faces in the grid   
!==============================================================================!

  ! Store the number of faces for the grid
  grid % n_faces = nf

  ! Number of nodes at each face (determines face's shape really)
  allocate(grid % faces_n_nodes(nf+ns));  grid % faces_n_nodes = 0

  ! Faces' nodes and neigboring cells
  allocate(grid % faces_n(4, nf+ns));  grid % faces_n = 0
  allocate(grid % faces_c(2, nf+ns));  grid % faces_c = 0

  ! Face surface areas (si) and distances between cells (di)
  allocate(grid % sx(nf+ns));  grid % sx = 0.0
  allocate(grid % sy(nf+ns));  grid % sy = 0.0
  allocate(grid % sz(nf+ns));  grid % sz = 0.0

  allocate(grid % dx(nf+ns));  grid % dx = 0.0
  allocate(grid % dy(nf+ns));  grid % dy = 0.0
  allocate(grid % dz(nf+ns));  grid % dz = 0.0

  ! Face center coordinates
  allocate(grid % xf(nf+ns));  grid % xf = 0.0
  allocate(grid % yf(nf+ns));  grid % yf = 0.0
  allocate(grid % zf(nf+ns));  grid % zf = 0.0

  ! Weight factors
  allocate(grid % f(nf+ns));   grid % f = 0.0

  end subroutine
