!==============================================================================!
  subroutine Allocate_Faces(Grid, nf, ns)
!------------------------------------------------------------------------------!
!>  Allocates memory for face-based data (arrays and matrices), for geometrical
!>  (xf, yf, zf, dx, dy, dz, ...) and connectivity data (faces_c, faces_n, ...).
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type)    :: Grid  !! grid under consideration
  integer, intent(in) :: nf    !! number of faces in the grid
  integer, intent(in) :: ns    !! number of shadow faces in the grid
!==============================================================================!

  ! Store the number of faces for the Grid
  Grid % n_faces   = nf
  Grid % n_shadows = ns

  ! Number of nodes at each face (determines face's shape really)
  allocate(Grid % faces_n_nodes(nf+ns));  Grid % faces_n_nodes(:) = 0

  ! Faces' nodes, neigboring cells and shadows
  allocate(Grid % faces_n(3, nf+ns));  Grid % faces_n(:,:) = 0
  allocate(Grid % faces_c(2, nf+ns));  Grid % faces_c(:,:) = 0
  allocate(Grid % faces_s(   nf+ns));  Grid % faces_s  (:) = 0

  ! Face surface areas (si), total surface (s)
  ! and distances between cells (di)
  allocate(Grid % sx(nf+ns));  Grid % sx(:) = 0.0
  allocate(Grid % sy(nf+ns));  Grid % sy(:) = 0.0
  allocate(Grid % sz(nf+ns));  Grid % sz(:) = 0.0

  allocate(Grid % s (nf+ns));  Grid % s (:) = 0.0

  allocate(Grid % dx(nf+ns));  Grid % dx(:) = 0.0
  allocate(Grid % dy(nf+ns));  Grid % dy(:) = 0.0
  allocate(Grid % dz(nf+ns));  Grid % dz(:) = 0.0
  allocate(Grid % d (nf+ns));  Grid % d (:) = 0.0

  ! Face center coordinates
  allocate(Grid % xf(nf+ns));  Grid % xf(:) = 0.0
  allocate(Grid % yf(nf+ns));  Grid % yf(:) = 0.0
  allocate(Grid % zf(nf+ns));  Grid % zf(:) = 0.0

  ! Vectors connecting face center with face cell centers connection
  allocate(Grid % rx(nf+ns));  Grid % rx(:) = 0.0
  allocate(Grid % ry(nf+ns));  Grid % ry(:) = 0.0
  allocate(Grid % rz(nf+ns));  Grid % rz(:) = 0.0

  ! Fractional cell volumes around faces
  allocate(Grid % dv1(nf+ns));  Grid % dv1(:) = 0.0
  allocate(Grid % dv2(nf+ns));  Grid % dv2(:) = 0.0

  ! Weight factors
  allocate(Grid % f (nf+ns));   Grid % f (:) = 0.0
  allocate(Grid % fw(nf+ns));   Grid % fw(:) = 0.0

  ! Allocate thread i.d.
  allocate(Grid % Omp % face_thread(nf+ns));  Grid % Omp % face_thread(:) = 0

  ! Allocate new and old numbers (this is so often used, maybe is better here)
  if(PROGRAM_NAME .ne. "Process") then
    allocate(Grid % new_f(nf+ns));  Grid % new_f(:) = 0
    allocate(Grid % old_f(nf+ns));  Grid % old_f(:) = 0
  end if

  end subroutine
