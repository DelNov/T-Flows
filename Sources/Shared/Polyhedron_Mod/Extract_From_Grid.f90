!==============================================================================!
  subroutine Extract_From_Grid(Pol, Grid, cell, phi_n)
!------------------------------------------------------------------------------!
!>  The Extract_From_Grid subroutine in Polyhedron_Mod is designed to convert
!>  a T-Flows cell from a computational grid into a polyhedron compatible with
!>  the Isoap library.  The subroutine takes a Polyhedron_Type object, a
!>  Grid_Type object, a cell identifier, and an optional real array phi_n as
!>  inputs. It uses various internal procedures and local variables to transform
!>  the grid cell into a polyhedron compatible with Isoap's requirements.
!>  The debug mode can be enabled for additional output files.
!------------------------------------------------------------------------------!
!   Functionaity                                                               !
!                                                                              !
!   * Allocates memory for the polyhedron if not already done.                 !
!   * Initializes the polyhedron by resetting its properties.                  !
!   * Extracts the number of nodes and faces from the specified cell.          !
!   * Copies node coordinates and optional nodal phi values to the             !
!     polyhedron.                                                              !
!   * Determines local node indices for each polyhedron cell, ensuring         !
!     correct node order for Isoap processing.                                 !
!   * Optionally, plots the extracted cell for verification.                   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Polyhedron_Type),  intent(out) :: Pol       !! parent class
  type(Grid_Type), target, intent(in)  :: Grid      !! grid object
  integer,                 intent(in)  :: cell      !! cell number
  real,          optional, intent(in)  :: phi_n(:)  !! values at cell's nodes
!------------------------------[Local parameters]------------------------------!
  logical, parameter :: DEBUG = .false.  ! if true, a lot of files are created
!-----------------------------------[Locals]-----------------------------------!
  integer, pointer, contiguous :: glo(:)
  integer, pointer, contiguous :: local_node(:)         ! local to polyhedron
  integer                      :: local_face_nodes(MAX_ISOAP_VERTS)
  integer                      :: i_nod, i_fac, l_nod
  integer                      :: s, n, faces_n_nodes
  real                         :: dx, dy, dz, sx, sy, sz
!==============================================================================!

  ! Take alias for global cell numbers
  glo => Grid % Comm % cell_glo

  call Work % Connect_Int_Node(local_node)  ! this also sets it to zero

  !----------------------------------------------------------------------!
  !   If not done yet, allocate memory for polyhedron and iso-polygons   !
  !----------------------------------------------------------------------!
  if(.not. Pol % allocated) then
    call Pol % Allocate_Polyhedron(MAX_ISOAP_FACES, MAX_ISOAP_VERTS)
  end if

  !-------------------------------!
  !   (Re)initialize Polyhedron   !
  !-------------------------------!
  Pol % n_nodes            = 0
  Pol % n_faces            = 0
  Pol % faces_n_nodes(:)   = 0
  Pol % faces_n      (:,:) = 0
  Pol % nodes_xyz    (:,:) = 0.0
  Pol % phi          (:)   = 0.0
  Pol % phi_iso            = 0.5
  Pol % global_node  (:)   = 0

  !---------------------------------------!
  !   Extract number of nodes and faces   !
  !---------------------------------------!
  Pol % n_faces = Grid % cells_n_faces(cell)
  Pol % n_nodes = abs(Grid % cells_n_nodes(cell))  ! < 0 for polyhedral

  !----------------------------------------------------------------------!
  !   Find local node indices and copy their coordinates and nodal phi   !
  !----------------------------------------------------------------------!
  l_nod = 0                                      ! initialize local node count
  do i_nod = 1, abs(Grid % cells_n_nodes(cell))  ! local (to cell) node number

    ! Form local node indices
    n = Grid % cells_n(i_nod, cell)         ! global (to grid) node number
    l_nod = l_nod + 1                       ! increase local node count
    local_node(n) = l_nod                   ! store local node number

    ! Copy node coordinates to polyhedron
    Pol % nodes_xyz(l_nod,1) = Grid % xn(n)
    Pol % nodes_xyz(l_nod,2) = Grid % yn(n)
    Pol % nodes_xyz(l_nod,3) = Grid % zn(n)

    ! Since you are here, copy the nodal phi values and global node numbers
    if(present(phi_n)) then
      Pol % phi(l_nod) = phi_n(n)
    end if
    Pol % global_node(l_nod) = n
  end do

  !------------------------------------------------------!
  !   Find local node indices for each polyhedron cell   !
  !------------------------------------------------------!
  do i_fac = 1, Grid % cells_n_faces(cell)    ! local (to cell) face number
    s  = Grid % cells_f(i_fac, cell)          ! global (to grid) face number
    if(Grid % faces_s(s) .ne. 0) then         ! if face has a shadow ...
      if(Grid % faces_c(1,s) .ne. cell) then  ! ... and cell is not its c1 ...
        s = Grid % faces_s(s)                 ! ... take the shadow instead
      end if
    end if

    Assert(Grid % Is_Face_In_Cell(s, cell))

    ! Fetch the local node numbers for each face
    faces_n_nodes = Grid % faces_n_nodes(s)  ! number of nodes in this face
    do i_nod = 1, faces_n_nodes              ! local (to face) node number
      n     = Grid % faces_n(i_nod, s)       ! global (to grid) node number
      l_nod = local_node(n)                  ! local (to polyhedron) node number
      local_face_nodes(i_nod) = l_nod        ! store in array of local nodes
    end do

    call Grid % Faces_Surface(s, sx, sy, sz)
    dx = Grid % xf(s) - Grid % xc(cell)
    dy = Grid % yf(s) - Grid % yc(cell)
    dz = Grid % zf(s) - Grid % zc(cell)

    ! Faces' nodes might be in the wrong order, correct here if needed.
    ! (If the face is oriented outwards to current cell,the nodes have to
    ! sorted in reverse order because ISOAP requies faces to point inwards.
    ! But then again, this seems to make the extracted iso-surface to point
    ! from VOF = 1 towards VOF = 0, which is not in line with T-Flows.)
    if(sx*dx + sy*dy + sz*dz < 0) then
      call Sort % Reverse_Order_Int(local_face_nodes(1:faces_n_nodes))
    end if

    ! Store everything in polyhedron
    Pol % faces_n_nodes(i_fac) = faces_n_nodes
    Pol % faces_n      (i_fac, 1:faces_n_nodes)  &
            = local_face_nodes(1:faces_n_nodes)
  end do

  ! Plot extracted cell first, in case things go wrong
  if(DEBUG) then
    ! This will create "geo-000XXX.vtk"
    call Pol % Plot_Polyhedron_Vtk("geo", glo(cell))
  end if

  call Work % Disconnect_Int_Node(local_node)

  end subroutine
