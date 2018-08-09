!==============================================================================!
  subroutine Grid_Mod_Allocate_Faces(grid, nf)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: nf    ! number of faces in the grid   
!==============================================================================!

  ! Store the number of faces for the grid
  grid % n_faces = nf

  ! Number of nodes at each face (determines face's shape really)
  allocate(grid % faces_n_nodes(nf));  grid % faces_n_nodes = 0

  ! Faces' nodes and neigboring cells
  allocate(grid % faces_n(4, nf));  grid % faces_n = 0
  allocate(grid % faces_c(2, nf));  grid % faces_c = 0

  ! Face surface areas (si), total surface (s) 
  ! and distances between cells (di)
  allocate(grid % sx(nf));  grid % sx = 0.0
  allocate(grid % sy(nf));  grid % sy = 0.0
  allocate(grid % sz(nf));  grid % sz = 0.0

  allocate(grid % s (nf));  grid % s  = 0.0

  allocate(grid % dx(nf));  grid % dx = 0.0
  allocate(grid % dy(nf));  grid % dy = 0.0
  allocate(grid % dz(nf));  grid % dz = 0.0

  ! Face center coordinates
  allocate(grid % xf(nf));  grid % xf = 0.0
  allocate(grid % yf(nf));  grid % yf = 0.0
  allocate(grid % zf(nf));  grid % zf = 0.0

  ! Weight factors
  allocate(grid % f(nf));   grid % f = 0.0

  end subroutine
