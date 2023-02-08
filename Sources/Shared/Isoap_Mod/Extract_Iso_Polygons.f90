!==============================================================================!
  subroutine Extract_Iso_Polygons(Isoap, Grid, cell, phi_n)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Isoap_Type)   :: Isoap
  type(Grid_Type)     :: Grid
  integer, intent(in) :: cell
  real                :: phi_n(:)
!------------------------------[Local parameters]------------------------------!
  logical, parameter :: DEBUG = .false.  ! if true, a lot of files are created
!-----------------------------------[Locals]-----------------------------------!
  integer, pointer, contiguous :: glo(:)
  integer                      :: local_face_nodes(MAX_ISOAP_VERTS)
  integer                      :: i_nod, i_fac, i_ver, i_iso, l_nod
  integer                      :: s, n, faces_n_nodes
!==============================================================================!

  ! Take alias for global cell numbers
  glo => Grid % Comm % cell_glo

  !-------------------------------------------------------------------------!
  !   On the first visit, allocate memory for polyhedron and iso-polygons   !
  !-------------------------------------------------------------------------!
  if(.not. Iso_Polygons % allocated) then
    call Iso_Polygons % Allocate_Iso_Polygons(MAX_ISOAP_FACES, MAX_ISOAP_VERTS)
  end if

  !------------------------!
  !   Extract polyhedron   !
  !------------------------!
  call Polyhedron % Extract_From_Grid(Grid, cell, phi_n)

  if(DEBUG) then
    ! This will create "extract-000XXX.vtk"
    call Polyhedron % Plot_Polyhedron_Vtk("extract", glo(cell))
  end if

  !---------------------------------!
  !   (Re)initialize Iso_Polygons   !
  !---------------------------------!
  Iso_Polygons % n_polys            = 0
  Iso_Polygons % polys_v      (:,:) = 0
  Iso_Polygons % face_index   (:)   = 0
  Iso_Polygons % polys_n_verts(:)   = 0
  Iso_Polygons % verts_xyz    (:,:) = 0.0
  Iso_Polygons % b_node_1     (:)   = 0
  Iso_Polygons % b_node_2     (:)   = 0

  !------------------------------!
  !   Call the Isoap algorithm   !
  !------------------------------!
  call Isoap % Main_Isoap(Polyhedron, Iso_Polygons)

  ! Plot extracted polygons
  if(DEBUG) then
     ! This will create "iso-000XXX.vtk"
     call Iso_Polygons % Plot_Iso_Polygons_Vtk(glo(cell))
  end if

  end subroutine
