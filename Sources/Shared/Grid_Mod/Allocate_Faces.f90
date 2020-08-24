!==============================================================================!
  subroutine Grid_Mod_Allocate_Faces(grid, nf, ns)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: nf    ! number of faces in the grid   
  integer         :: ns    ! number of shadow faces in the grid   
!-----------------------------------[Locals]-----------------------------------!
  integer         :: s
!==============================================================================!

  ! Store the number of faces for the grid
  grid % n_faces   = nf
  grid % n_shadows = ns

  ! Number of nodes at each face (determines face's shape really)
  allocate(grid % faces_n_nodes(nf+ns));  grid % faces_n_nodes = 0

  ! Faces' nodes, neigboring cells and shadows
  allocate(grid % faces_n(4,nf+ns));  grid % faces_n = 0
  allocate(grid % faces_c(2,nf+ns));  grid % faces_c = 0
  allocate(grid % faces_s(  nf+ns));  grid % faces_s = 0

  ! Face surface areas (si), total surface (s) 
  ! and distances between cells (di)
  allocate(grid % sx(nf+ns));  grid % sx = 0.0
  allocate(grid % sy(nf+ns));  grid % sy = 0.0
  allocate(grid % sz(nf+ns));  grid % sz = 0.0

  allocate(grid % s (nf+ns));  grid % s  = 0.0

  allocate(grid % dx(nf+ns));  grid % dx = 0.0
  allocate(grid % dy(nf+ns));  grid % dy = 0.0
  allocate(grid % dz(nf+ns));  grid % dz = 0.0
  allocate(grid % d(nf+ns));   grid % d = 0.0

  ! Face center coordinates
  allocate(grid % xf(nf+ns));  grid % xf = 0.0
  allocate(grid % yf(nf+ns));  grid % yf = 0.0
  allocate(grid % zf(nf+ns));  grid % zf = 0.0

  ! Vectors connecting face center with face cell centers connection
  allocate(grid % xr(nf+ns));  grid % xr = 0.0
  allocate(grid % yr(nf+ns));  grid % yr = 0.0
  allocate(grid % zr(nf+ns));  grid % zr = 0.0

  ! Weight factors
  allocate(grid % f (nf+ns));   grid % f  = 0.0
  allocate(grid % fw(nf+ns));   grid % fw = 0.0

  ! Allocate processor i.d.
  allocate(grid % comm % face_glo(nf+ns));  grid % comm % face_glo(:) = 0
  do s = 1, nf+ns
    grid % comm % face_glo(s) = s
  end do

  end subroutine
