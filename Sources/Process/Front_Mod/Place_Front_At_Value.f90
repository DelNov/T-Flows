!==============================================================================!
  subroutine Place_Front_At_Value(Front, sharp, verbose)
!------------------------------------------------------------------------------!
!>  This subroutine is essential for defining a front in simulations,
!>  particularly for phase interfaces in VOF methods. It constructs the front
!>  where a specified field variable, typically the VOF function, reaches a
!>  value of 0.5, indicating the interface between two phases.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Initializing the front structure for fresh computation.                  !
!   * Interpolating cell values to nodes.                                      !
!   * Iterating through all cells to locate surface vertices.                  !
!   * Employing the Isoap algorithm to extract polygons representing the       !
!     interface.                                                               !1
!   * Storing results from Isoap into T-Flows' front object.                   !
!   * Compressing front vertices for efficient representation and less memory  !
!     usage.                                                                   !
!   * Establishing connectivity among elements and sides of the front.         !
!   * Calculating geometric quantities like centroids and normals of elements. !
!   * Marking cells and intersecting faces with elements for phase transfer    !
!     calculations.                                                            !
!   * Optionally saving debug information for visualization.                   !
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Front_Type), target :: Front    !! parent class
  type(Var_Type),    target :: sharp    !! variable for front placement,
                                        !! typically a sharp VOF function
  logical                   :: verbose  !! controls the output verbosity
!------------------------------[Local parameters]------------------------------!
  logical, parameter :: DEBUG = .false.
!-----------------------------------[Locals]-----------------------------------!
  type(Polyhedron_Type),   pointer :: Pol
  type(Iso_Polygons_Type), pointer :: Iso
  type(Grid_Type),         pointer :: Grid
  type(Field_Type),        pointer :: Flow
  type(Vert_Type),         pointer :: Vert(:)
  type(Elem_Type),         pointer :: Elem(:)
  integer,                 pointer :: nv, ne
  integer                          :: c, i, nb, nc, n, nn
  integer                          :: v, n_verts_in_buffers
  integer                          :: i_nod, i_ver, i_iso
  real                             :: phi_cell_min, phi_cell_max
  real, contiguous,        pointer :: phi_n(:)
!==============================================================================!

  call Profiler % Start('Creating_Front_From_Vof_Function')

  call Work % Connect_Real_Node(phi_n)

  ! Take aliases
  Pol  => Polyhedron
  Iso  => Iso_Polygons
  Grid => Front % pnt_grid
  Flow => Front % pnt_flow
  nv   => Front % n_verts
  ne   => Front % n_elems
  Vert => Front % Vert
  Elem => Front % Elem
  nb   =  Grid % n_bnd_cells
  nc   =  Grid % n_cells
  nn   =  Grid % n_nodes

  call Front % Initialize_Front()
  call Flow % Interpolate_Cells_To_Nodes(sharp % n, phi_n(1:nn))

  if(DEBUG) then
    call Grid % Save_Debug_Vtu('phi_c',                &
                               scalar_cell=sharp % n,  &
                               scalar_name='phi_c')
    call Grid % Save_Debug_Vtu('phi_n',              &
                               scalar_node=phi_n,    &
                               scalar_name='phi_n')
  end if

  !-----------------------------------------------!
  !   This is a bit ad-hoc - if some points are   !
  !    "exactly" 0.5, change them a little bit    !
  !-----------------------------------------------!
  do c = 1, nc
    do i_nod = 1, abs(Grid % cells_n_nodes(c))
      n = Grid % cells_n(i_nod, c)
      if(Math % Approx_Real(phi_n(n), 0.5)) then
        if(sharp % n(c) < sharp % o(c)) then
          phi_n(n) = phi_n(n) - MILI
        else
          phi_n(n) = phi_n(n) + MILI
        end if
      end if
    end do
  end do

  ! Find vertices in each cell
  nv = 0
  ne = 0

  !------------------------------------------------------------!
  !                                                            !
  !   Browse through all cells in search of surface vertices   !
  !                                                            !
  !------------------------------------------------------------!
  do c = Cells_In_Domain()
    phi_cell_min =  HUGE
    phi_cell_max = -HUGE
    do i_nod = 1, abs(Grid % cells_n_nodes(c))
      n = Grid % cells_n(i_nod, c)
      phi_cell_min = min(phi_cell_min, phi_n(n))
      phi_cell_max = max(phi_cell_max, phi_n(n))
    end do
    if(phi_cell_min <= 0.5 .and. phi_cell_max >= 0.5) then

      !------------------------------!
      !   Call the Isoap algorithm   !
      !------------------------------!
      call Isoap % Extract_Iso_Polygons(Grid, c, phi_n)

      !------------------------------------------------------!
      !   Store results from Isoap into T-Flows' variables   !
      !------------------------------------------------------!
      do i_iso = 1, Iso % n_polys

        ! Increase element's count
        ne = ne + 1

        Elem(ne) % cell = c

        Elem(ne) % nv = Iso % polys_n_verts(i_iso)

        ! Store elements vertices
        do i_ver = 1, Iso % polys_n_verts(i_iso)
          Elem(ne) % v(i_ver) = nv + Iso % polys_v(i_iso, i_ver)
        end do

        ! Fetch coordinates and bounding nodes of new iso-polygon
        do i_ver = 1, Iso % polys_n_verts(i_iso)
          i = Iso % polys_v(i_iso, i_ver)

          Vert(nv+i) % x_n = Iso % verts_xyz(i, 1)
          Vert(nv+i) % y_n = Iso % verts_xyz(i, 2)
          Vert(nv+i) % z_n = Iso % verts_xyz(i, 3)

          Front % b_node_1(nv+i) = Pol % global_node(Iso % b_node_1(i))
          Front % b_node_2(nv+i) = Pol % global_node(Iso % b_node_2(i))
        end do
      end do    ! through new polygons

      ! Update the number of new vertices
      do i_iso = 1, Iso % n_polys
        nv = nv + Iso % polys_n_verts(i_iso)
      end do

    end if
  end do

  !-----------------------!
  !                       !
  !   Compress vertices   !
  !                       !
  !-----------------------!
  call Front % Compress_Front_Vertices(verbose)

  !---------------------------------!
  !                                 !
  !   Find and check connectivity   !
  !                                 !
  !---------------------------------!
  call Front % Find_Sides(verbose)    ! Surf calls the same here
  call Front % Find_Front_Elements()  ! Surf calls Find_Surf_Elements
  call Front % Check_Elements()       ! Surf calls the same here

  !-----------------------------------------------!
  !   It used to find the nearest cell and node   !
  !   (But I deleted, I hope it was not needed)   !
  !-----------------------------------------------!
  n_verts_in_buffers = 0
  do v = 1, nv
    call Front % Vert(v) % Find_Nearest_Cell(n_verts_in_buffers, locally=.true.)
    call Front % Vert(v) % Find_Nearest_Node()
  end do

  !--------------------------------------!
  !                                      !
  !   Calculate geometrical quantities   !
  !                                      !
  !--------------------------------------!
  call Front % Find_Vertex_Elements()
  call Front % Calculate_Element_Centroids()
  call Front % Calculate_Element_Normals()

  !-----------------------------------------------------------------!
  !                                                                 !
  !    Find cells at elements and intersect faces with elements     !
  !   (This is needed only for mass transfer problems with front)   !
  !                                                                 !
  !-----------------------------------------------------------------!
  call Front % Mark_Cells_And_Faces(sharp)

  if(DEBUG) then
    call Front % Save_Debug_Front_Vtu(0)  ! 0 is for time step
  end if

  call Work % Disconnect_Real_Node(phi_n)

  call Profiler % Stop('Creating_Front_From_Vof_Function')

  return

  end subroutine
