!==============================================================================!
  subroutine Extract_Iso_Polygons(Isoap, Grid, cell, phi_n)
!------------------------------------------------------------------------------!
!>  The Extract_Iso_Polygons subroutine is pivotal in T-Flows' interaction with
!>  Isoap. It extracts a polyhedron from the computational grid using the
!>  Polyhedron % Extract_From_Grid function and stores it in the singleton
!>  object Polyhedron.  After that, it processes it through Isoap's algorithm
!>  via Isoap % Main_Isoap, resulting in iso-polygons stored in singleton
!>  object Iso_Polygones.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Isoap_Type),       intent(out) :: Isoap     !! parent class
  type(Grid_Type), target, intent(in)  :: Grid      !! computational grid
  integer,                 intent(in)  :: cell      !! cell being analyze
  real,          optional, intent(in)  :: phi_n(:)  !! variable interpolated
                                                    !! at the nodes
!------------------------------[Local parameters]------------------------------!
  logical, parameter :: DEBUG = .false.  ! if true, a lot of files are created
!-----------------------------------[Locals]-----------------------------------!
  integer, pointer, contiguous :: glo(:)
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
    ! This will create "extract-cell-000XXX.vtk"
    call Polyhedron % Plot_Polyhedron_Vtk("extract-cell", glo(cell))
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
     ! This will create "extract-iso-000XXX.vtk"
     call Iso_Polygons % Plot_Iso_Polygons_Vtk("extract-iso", glo(cell))
  end if

  end subroutine
