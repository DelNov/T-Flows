!==============================================================================!
  subroutine Extract_From_Grid(Polyhedron, Grid, cell, phi_n)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Polyhedron_Type) :: Polyhedron
  type(Grid_Type)        :: Grid
  integer, intent(in)    :: cell
  real,    optional      :: phi_n(:)
!-----------------------------------[Locals]-----------------------------------!
  logical, save                :: first_visit = .true.
  integer, contiguous, pointer :: local_node(:)         ! local to polyhedron
  integer                      :: local_face_nodes(MAX_ISOAP_VERTS)
  integer                      :: i_nod, i_fac, i_ver, i_iso, l_nod
  integer                      :: s, n, faces_n_nodes
!==============================================================================!

  call Work % Connect_Int_Node(local_node)  ! this also sets it to zero

  !-------------------------------------------------------------------------!
  !   On the first visit, allocate memory for polyhedron and iso-polygons   !
  !-------------------------------------------------------------------------!
  if(first_visit) then
    call Polyhedron % Allocate_Polyhedron(MAX_ISOAP_FACES, MAX_ISOAP_VERTS)
    first_visit = .false.
  end if

  !------------------------------------------------!
  !   (Re)initialize Polyhedron and Iso_Polygons   !
  !------------------------------------------------!
  Polyhedron % n_nodes            = 0
  Polyhedron % n_faces            = 0
  Polyhedron % faces_n_nodes(:)   = 0
  Polyhedron % faces_n      (:,:) = 0
  Polyhedron % nodes_xyz    (:,:) = 0.0
  Polyhedron % phi          (:)   = 0.0
  Polyhedron % phiiso             = 0.5
  Polyhedron % global_node  (:)   = 0

  !---------------------------------------!
  !   Extract number of nodes and faces   !
  !---------------------------------------!
  Polyhedron % n_faces = Grid % cells_n_faces(cell)
  Polyhedron % n_nodes = abs(Grid % cells_n_nodes(cell))  ! < 0 for polyhedral

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
    Polyhedron % nodes_xyz(l_nod,1) = Grid % xn(n)
    Polyhedron % nodes_xyz(l_nod,2) = Grid % yn(n)
    Polyhedron % nodes_xyz(l_nod,3) = Grid % zn(n)

    ! Since you are here, copy the nodal phi values and global node numbers
    if(present(phi_n)) then
      Polyhedron % phi(l_nod) = phi_n(n)
    end if
    Polyhedron % global_node(l_nod) = n
  end do

  !------------------------------------------------------!
  !   Find local node indices for each polyhedron cell   !
  !------------------------------------------------------!
  do i_fac = 1, Grid % cells_n_faces(cell)  ! local (to cell) face number
    s = Grid % cells_f(i_fac, cell)         ! global (to grid) face number

    ! Fetch the local node numbers for each face
    faces_n_nodes = Grid % faces_n_nodes(s)  ! number of nodes in this face
    do i_nod = 1, faces_n_nodes              ! local (to face) node number
      n     = Grid % faces_n(i_nod, s)       ! global (to grid) node number
      l_nod = local_node(n)                  ! local (to polyhedron) node number
      local_face_nodes(i_nod) = l_nod        ! store in array of local nodes
    end do

    ! They might be in the wrong order, correct if needed
    ! (If the face is oriented inwards to current cell,
    ! the nodes have to sorted in reverse order.)
    if(Grid % cells_f_orient(i_fac, cell) .eq. OUTWARDS) then
      call Sort % Reverse_Order_Int(local_face_nodes(1:faces_n_nodes))
    end if

    ! Store everything in polyhedron
    Polyhedron % faces_n_nodes(i_fac) = faces_n_nodes
    Polyhedron % faces_n      (i_fac, 1:faces_n_nodes)  &
                   = local_face_nodes(1:faces_n_nodes)
  end do

  ! Plot extracted cell first, in case things go wrong
  ! call Polyhedron % Plot_Polyhedron_Vtk(cell)

  call Work % Disconnect_Int_Node(local_node)

  end subroutine
