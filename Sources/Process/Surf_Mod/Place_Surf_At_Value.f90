!==============================================================================!
  subroutine Place_Surf_At_Value(Surf, sharp, smooth, verbose)
!------------------------------------------------------------------------------!
!>  This subroutine is crucial for defining the position of the surface within
!>  a computational domain based on a specified value of a scalar field. It
!>  operates by identifying intersection points within the computational grid
!>  where the scalar field value transitions through the specified threshold,
!>  facilitating the creation of a detailed and accurate representation of the
!>  fluid interface.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Initialization and connection: Connects to the real node of the scalar   !
!     field and initializes the surface.                                       !
!   * Alias setup and interpolation: Sets up aliases for grid, flow, vertices, !
!     and elements, and interpolates cell values to nodes.                     !
!   * Gradient calculation: Computes gradients of the scalar field for more    !
!     accurate surface placement.                                              !
!   * Vertex identification and placement:                                     !
!     - Iterates through grid cells, identifying vertices where the scalar     !
!       field intersects the specified value.                                  !
!     - Calculates the precise position of these vertices based on linear      !
!       interpolation between grid nodes.                                      !
!   * Element formation:                                                       !
!     * Determines the number of vertices intersecting in each cell and forms  !
!       elements based on these vertices.                                      !
!     * Handles different cases (3-point, 4-point, etc.) using specific        !
!       subroutines for each scenario.                                         !
!     * Ensures that elements formed are valid and consistent.                 !
!   * Parallel processing considerations: In parallel runs, distributes the    !
!     mesh to ensure consistency across processors using Distribute_Mesh.      !
!   * Vertex compression: Compresses surface vertices to eliminate duplicates  !
!     and optimize the mesh.                                                   !
!   * Connectivity and verification: Establishes sides and elements, checks    !
!     for consistency and accuracy in the mesh.                                !
!   * Nearest cell and node identification: For each vertex, finds the nearest !
!     cell and node in the computational grid.                                 !
!   * Smooth function value and gradient storage: Distributes the smooth       !
!     function value and its gradients to vertices.                            !
!   * Geometrical calculations: Computes centroids and normals for elements,   !
!     enhancing the geometric accuracy of the mesh.                            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Surf_Type),  target :: Surf     !! parent class
  type(Var_Type),    target :: sharp    !! sharp scalar function (mostly VOF)
  type(Var_Type),    target :: smooth   !! smoothed scalar function
  logical                   :: verbose  !! controls output verbosity
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: Grid
  type(Field_Type),  pointer :: Flow
  type(Vert_Type),   pointer :: Vert(:)
  type(Elem_Type),   pointer :: Elem(:)
  integer,           pointer :: nv, ne
  integer, allocatable       :: n_cells_v(:)
  integer                    :: c, j, n1, n2, nb, nc, nn
  integer                    :: v, n_vert, n_verts_in_buffers
  integer                    :: en(12,2)  ! edge numbering
  real                       :: phi1, phi2, xn1, yn1, zn1, xn2, yn2, zn2, w1, w2
  real                       :: surf_v(3)
  real, contiguous,  pointer :: phi_n(:)
!------------------------------------------------------------------------------!
  include 'Surf_Mod/Edge_Numbering.h90'
!==============================================================================!

  call Work % Connect_Real_Node(phi_n)

  call Profiler % Start('Creating_Surface_From_Vof_Function')

  ! Take aliases
  Grid => Surf % pnt_grid
  Flow => Surf % pnt_flow
  nv   => Surf % n_verts
  ne   => Surf % n_elems
  Vert => Surf % Vert
  Elem => Surf % Elem
  nb   =  Grid % n_bnd_cells
  nc   =  Grid % n_cells
  nn   =  Grid % n_nodes

  call Surf % Initialize_Surf()
  call Flow % Interpolate_Cells_To_Nodes(sharp % n, phi_n(1:nn))

  !---------------------------------------------!
  !   Take gradients from the smooth function   !
  !   These might be already computed - check   !
  !---------------------------------------------!
  call Flow % Grad_Variable(smooth)

  allocate(n_cells_v(Grid % n_cells))
  n_cells_v(:) = 0

  ! Find vertices in each cell
  nv = 0
  ne = 0

  !------------------------------------------------------------!
  !                                                            !
  !   Browse through all cells in search of surface vertices   !
  !                                                            !
  !------------------------------------------------------------!
  do c = 1, Grid % n_cells - Grid % Comm % n_buff_cells

    n_vert = 0

    ! Fetch the edges for this cell
    if(Grid % cells_n_nodes(c) .eq. 4) en = edg_tet
    if(Grid % cells_n_nodes(c) .eq. 5) en = edg_pyr
    if(Grid % cells_n_nodes(c) .eq. 6) en = edg_wed
    if(Grid % cells_n_nodes(c) .eq. 8) en = edg_hex

    !------------------------------------------------------!
    !   Browse through edges to find intersection points   !
    !------------------------------------------------------!
    do j = 1, 12  ! max number of edges
      n1 = Grid % cells_n( en(j,1), c )
      n2 = Grid % cells_n( en(j,2), c )

      phi1 = phi_n(n1)
      phi2 = phi_n(n2)

      ! There is a vertex between these two edges
      if( ((phi2 - 0.5) * (0.5 - phi1)) >= MICRO ) then
        n_cells_v(c) = n_cells_v(c) + 1

        nv = nv + 1

        n_vert = n_vert + 1

        xn1 = Grid % xn(n1)
        yn1 = Grid % yn(n1)
        zn1 = Grid % zn(n1)

        xn2 = Grid % xn(n2)
        yn2 = Grid % yn(n2)
        zn2 = Grid % zn(n2)

        w1 = abs(phi2 - 0.5) / abs(phi2 - phi1)
        w2 = 1.0 - w1

        ! All vertices have to be stored
        Vert(nv) % x_n = xn1*w1 + xn2*w2
        Vert(nv) % y_n = yn1*w1 + yn2*w2
        Vert(nv) % z_n = zn1*w1 + zn2*w2

      end if

    end do  ! through edges

    !---------------------------------!
    !   Some points have been found   !
    !---------------------------------!
    if(n_vert > 0) then

      ! Surface vector
      surf_v(1) = smooth % x(c)
      surf_v(2) = smooth % y(c)
      surf_v(3) = smooth % z(c)

      ! If valid elements were formed (last argument: enforce_triangles)
      if(n_vert .eq. 3) call Surf % Handle_3_Points(surf_v)
      if(n_vert .eq. 4) call Surf % Handle_4_Points(surf_v, .true.)
      if(n_vert .eq. 5) call Surf % Handle_5_Points(surf_v, .true.)
      if(n_vert .eq. 6) call Surf % Handle_6_Points(surf_v, .true.)
      if(n_vert .eq. 7) then
        print *, '# ERROR: seven vertices in an intersection!'
        stop
      end if
    end if

  end do

  !-----------------------------------------------------------------------!
  !                                                                       !
  !   At this point, each processor has its own vertices and some of      !
  !   them are, of course, duplicate.  Elements, on the other hand,       !
  !   should be unique since buffer cells are avoided in the loop above   !
  !                                                                       !
  !-----------------------------------------------------------------------!
  if(Parallel_Run()) then
    call Surf % Distribute_Mesh(verbose)
  end if

  !-----------------------!
  !                       !
  !   Compress vertices   !
  !                       !
  !-----------------------!
  call Surf % Compress_Surf_Vertices(verbose)

  !---------------------------------!
  !                                 !
  !   Find and check connectivity   !
  !                                 !
  !---------------------------------!
  call Surf % Find_Sides(verbose)   ! Front calls the same here
  call Surf % Find_Surf_Elements()  ! Front calls Find_Front_Elements
  call Surf % Check_Elements()      ! Front calls the same here

  !--------------------------------!
  !   Find nearest cell and node   !
  !--------------------------------!
  n_verts_in_buffers = 0
  do v = 1, nv
    call Surf % Vert(v) % Find_Nearest_Cell(n_verts_in_buffers)
    call Surf % Vert(v) % Find_Nearest_Node()
  end do

  !--------------------------------------------!
  !   Store the value of the smooth function   !
  !      and its gradients in the vertex       !
  !--------------------------------------------!
  call Surf % Distribute_Smooth(smooth)

  !--------------------------------------!
  !                                      !
  !   Calculate geometrical quantities   !
  !                                      !
  !--------------------------------------!
  call Surf % Find_Vertex_Elements()
  call Surf % Calculate_Element_Centroids()
  call Surf % Calculate_Element_Normals()

  call Profiler % Stop('Creating_Surface_From_Vof_Function')

  call Work % Disconnect_Real_Node(phi_n)

  return

  end subroutine
