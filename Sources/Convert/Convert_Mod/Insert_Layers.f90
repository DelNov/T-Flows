!==============================================================================!
  subroutine Insert_Layers(Convert, Grid,                              &
                                    n_layer_regions, layer_regions,    &
                                    n_layers,        layer_thickness,  &
                                    compress_domain)
!------------------------------------------------------------------------------!
!   At this point, boundary and inside faces have been found, but
!   geometrical quantities and final sorted ranges have not been formed yet.
!
!   Reliable:
!     Grid % n_nodes, n_cells, n_bnd_cells, n_faces, n_shadows
!     Grid % xn
!     Grid % yn
!     Grid % zn
!     Grid % cells_n_nodes, cells_n
!     Grid % faces_n_nodes, faces_n, faces_c
!     Grid % region % at_cell(c2), for boundary c2 < 0
!
!   Allocated, but not yet calculated/reliable:
!     Grid % sx, sy, sz, s
!     Grid % xf, yf, zf
!     Grid % dx, dy, dz, d
!     Grid % rx, ry, rz
!     Grid % f, fw
!     Grid % xc, yc, zc, vol
!
! Do not rely yet on:
!   Grid % cells_f, cells_n_faces, cells_c
!   boundary faces being in final region-sorted ranges
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Convert_Type)     :: Convert  !! parent class
  type(Grid_Type), target :: Grid     !! primal grid
  integer,    intent(in)  :: n_layer_regions
  integer                 :: layer_regions(n_layer_regions)
  integer,    intent(in)  :: n_layers
  real                    :: layer_thickness(n_layers)
  logical,    intent(in)  :: compress_domain
!------------------------------[Local parameters]------------------------------!
  logical, parameter :: DEBUG = .false.
!-----------------------------------[Locals]-----------------------------------!
  integer              :: i_reg, reg, s, c1, c2, ni, nj
  integer              :: cnt_c, cnt_f, cnt_n, cnt_e, cnt_ei, cnt_eb, run
  integer              :: n, nn, nn1, nn2
  integer              :: i_nod, j_nod, e, ee, n1, n2, s1, s2, sm, sb
  integer              :: old_n, old_f, old_c, old_bc, cnt_bnd_n
  real                 :: area, dot, eps, magn
  real                 :: nx, ny, nz
  real                 :: xf, yf, zf, xc, yc, zc, sx, sy, sz
  real                 :: xc1, yc1, zc1, xc2, yc2, zc2, dx, dy, dz
  real                 :: thickness, a(3,3), b(3), x(3)
  logical              :: invertible, inner, outer, warning
  integer              :: i_lay
  integer, allocatable :: node_to(:), face_to_face(:), face_to_cell(:)
  integer, allocatable :: edge_cnt(:), key(:)
  integer, allocatable :: edge_n1(:), edge_n2(:), edge_s1(:), edge_s2(:)
  integer, allocatable :: mark_nodes(:), mark_face(:), node_cnt(:)
  logical, allocatable :: fixed_nodes(:)
  real,    allocatable :: phi(:)
  real,    allocatable :: disp_x(:), disp_y(:), disp_z(:)  ! node displacements
  real,    allocatable :: node_nx(:), node_ny(:), node_nz(:)  ! node normals
  real,    allocatable :: a11(:), a12(:), a13(:)
  real,    allocatable :: a22(:), a23(:), a33(:)
  real,    allocatable :: bx(:),  by(:),  bz(:)
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Convert)
!==============================================================================!

  do i_lay = 1, n_layers

    if(compress_domain) then
      thickness = sum(layer_thickness(1:n_layers))
    else
      thickness = layer_thickness(i_lay)
      print '(a,es10.3)', " # Inserting a boundary layer with thickness ",  &
                            thickness
    end if

    !-----------------------------------------------------------------------!
    !                                                                       !
    !   Prepare boundary nodes and normals for one inserted layer           !
    !                                                                       !
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
    !   All boundary nodes are moved to the beginning of the node arrays.   !
    !   This keeps boundary-only work compact, while still allowing the     !
    !   selected regions to use only the subset of nodes which are really   !
    !   extruded in the current layer.                                      !
    !                                                                       !
    !   The selected boundary faces are then marked and their outward       !
    !   unit normals are accumulated at the extruded nodes.  For each such  !
    !   node a small least-squares system is formed.  Its solution gives    !
    !   the node displacement direction which has unit projection on all    !
    !   adjacent marked-face normals.                                       !
    !-----------------------------------------------------------------------!

    !----------------------------------------------!
    !   Sort nodes boundaries first because only   !
    !   these will be coppied to new layer later   !
    !----------------------------------------------!
    call Convert % Sort_Nodes_Boundary_First(Grid, cnt_bnd_n)

    !--------------------------------------------------------!
    !   Count nodes which will be copied into the new layer  !
    !    in order to know how much more memory to allocate.  !
    !--------------------------------------------------------!
    call Enlarge % Array_Int (mark_nodes, (/1, Grid % n_nodes/), val=0)

    cnt_n = 0  ! number of extruded nodes
    do i_reg = 1, n_layer_regions
      reg = layer_regions(i_reg)

      ! Browse through all faces
      do s = 1, Grid % n_faces
        c2 = Grid % faces_c(2,s)

        ! Check if face belings to the selected region ...
        if(c2 < 0) then
          if(Grid % region % at_cell(c2) .eq. reg) then

            ! ... and if so, mark the face's nodes
            do i_nod = 1, Grid % faces_n_nodes(s)
              n = Grid % faces_n(i_nod, s)
              if(mark_nodes(n) .eq. 0) then
                cnt_n = cnt_n + 1  ! increase the node count
                mark_nodes(n) = 1  ! just raise the flag
              end if
            end do
          end if
        end if
      end do  ! through faces
    end do    ! through regions

    print *, "# Old number of nodes:             ", Grid % n_nodes
    print *, "# New nodes in the boundary layer: ", cnt_n

    !----------------------------------------------!
    !   Allocate memory for local working arrays   !
    !----------------------------------------------!
    call Enlarge % Array_Int (mark_face, (/1, Grid % n_faces/), val=NO)

    call Enlarge % Array_Real(disp_x, (/1, Grid % n_nodes/), val=0.0)
    call Enlarge % Array_Real(disp_y, (/1, Grid % n_nodes/), val=0.0)
    call Enlarge % Array_Real(disp_z, (/1, Grid % n_nodes/), val=0.0)

    call Enlarge % Array_Real(a11, (/1, cnt_bnd_n/), val=0.0)
    call Enlarge % Array_Real(a12, (/1, cnt_bnd_n/), val=0.0)
    call Enlarge % Array_Real(a13, (/1, cnt_bnd_n/), val=0.0)
    call Enlarge % Array_Real(a22, (/1, cnt_bnd_n/), val=0.0)
    call Enlarge % Array_Real(a23, (/1, cnt_bnd_n/), val=0.0)
    call Enlarge % Array_Real(a33, (/1, cnt_bnd_n/), val=0.0)
    call Enlarge % Array_Real(bx,  (/1, cnt_bnd_n/), val=0.0)
    call Enlarge % Array_Real(by,  (/1, cnt_bnd_n/), val=0.0)
    call Enlarge % Array_Real(bz,  (/1, cnt_bnd_n/), val=0.0)

    mark_nodes(:) = 0

    !-----------------------------------------------------------------------!
    !                                                                       !
    !   Mark faces in extruding regions and assemble nodal normal systems   !
    !                                                                       !
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
    !   Selected boundary faces are marked only once, even if a region      !
    !   appears more than once in the input.  Their surface vectors are     !
    !   normalized and oriented away from the adjacent inner cell.  Each    !
    !   face normal then contributes one row to the least-squares system    !
    !   stored at every node of the face.                                   !
    !-----------------------------------------------------------------------!
    warning = .false.
    do i_reg = 1, n_layer_regions
      reg = layer_regions(i_reg)

      do s = 1, Grid % n_faces
        c1 = Grid % faces_c(1,s)
        c2 = Grid % faces_c(2,s)

        if(c2 < 0) then
          if(Grid % region % at_cell(c2) .eq. reg) then

            if(mark_face(s) .eq. NO) then
              mark_face(s) = YES

              ! Calculate face's surface normal
              call Grid % Faces_Surface(s, sx, sy, sz)

              area = sqrt(sx**2 + sy**2 + sz**2)
              Assert(area > TINY)

              nx = sx / area
              ny = sy / area
              nz = sz / area

              call Grid % Faces_Center(s, xf, yf, zf)

              call Grid % Cells_Center(c1, xc, yc, zc)

              dot = nx*(xf-xc) + ny*(yf-yc) + nz*(zf-zc)
              if(dot < 0.0) then
                warning = .true.
                nx = -nx
                ny = -ny
                nz = -nz
              end if

              do i_nod = 1, Grid % faces_n_nodes(s)
                n = Grid % faces_n(i_nod, s)
                Assert(n .gt. 0)
                Assert(n .le. cnt_bnd_n)

                mark_nodes(n) = mark_nodes(n) + 1

                a11(n) = a11(n) + nx * nx
                a12(n) = a12(n) + nx * ny
                a13(n) = a13(n) + nx * nz
                a22(n) = a22(n) + ny * ny
                a23(n) = a23(n) + ny * nz
                a33(n) = a33(n) + nz * nz

                bx(n) = bx(n) + nx
                by(n) = by(n) + ny
                bz(n) = bz(n) + nz
              end do
            end if
          end if
        end if
      end do
    end do

    if(warning) then
      call Message % Warning(60,                                             &
        "At some boundary faces, the dot product between surface vector "//  &
        "and connection between face center and inner cell center was "  //  &
        "negative.  Although not catastrophic, it shouldn't occur at "   //  &
        "this stage of grid conversion", file = __FILE__, line = __LINE__)
    end if

    eps = MILI

    do n = 1, Grid % n_nodes
      if(mark_nodes(n) .gt. 0) then

        a(1,1) = a11(n) + eps
        a(1,2) = a12(n)
        a(1,3) = a13(n)
        a(2,1) = a12(n)
        a(2,2) = a22(n) + eps
        a(2,3) = a23(n)
        a(3,1) = a13(n)
        a(3,2) = a23(n)
        a(3,3) = a33(n) + eps

        b(1) = bx(n)
        b(2) = by(n)
        b(3) = bz(n)
        x(:) = 0.0

        call Math % Gaussian_Elimination(3, a, b, x, invertible)

        if(invertible) then
          disp_x(n) = x(1)
          disp_y(n) = x(2)
          disp_z(n) = x(3)
        else
          disp_x(n) = bx(n) / mark_nodes(n)
          disp_y(n) = by(n) / mark_nodes(n)
          disp_z(n) = bz(n) / mark_nodes(n)
        end if

      end if  ! mark_nodes(n) .gt. 0
    end do    ! through nodes

    deallocate(a11)
    deallocate(a12)
    deallocate(a13)
    deallocate(a22)
    deallocate(a23)
    deallocate(a33)
    deallocate(bx)
    deallocate(by)
    deallocate(bz)

    !---------------------!
    !                     !
    !   Compress domain   !
    !                     !
    !---------------------!
    if(compress_domain) then

      ! Reserve memory for fixed_nodes, array which holds boundary conditions
      ! This is needed for Solve_Node_Potential and Smooth_Node_Displacements
      call Enlarge % Array_Log(fixed_nodes, (/1, Grid % n_nodes/), val=.false.)

      ! Find edges.
      ! This is needed for Solve_Node_Potential and Smooth_Node_Displacements
      call Grid % Find_Edges()

      !---------------------------------------!
      !   Solve for distance-like potential   !
      !---------------------------------------!
      call Enlarge % Array_Real(phi, (/1, Grid % n_nodes/), val=0.0)

      ! Set boundary conditions for potential
      do n = 1, Grid % n_nodes
        if(mark_nodes(n) .gt. 0) then
          fixed_nodes(n) = .true.
        end if
      end do    ! through nodes

      call Convert % Solve_Node_Potential(Grid, fixed_nodes, phi, 1.0e-6, 1000)

      call Grid % Save_Debug_Vtu(append="node-potential",  &
                                 scalar_node = phi,        &
                                 scalar_name = "node-potential")

      ! Reset fixed nodes (this probably not needed)
      fixed_nodes(:) = .false.

      !----------------------------------------------------------------!
      !   Form fixed_nodes arrays for distance like potential solver   !
      !----------------------------------------------------------------!

      ! Set boundary conditions for displacements
      do n = 1, Grid % n_nodes

        disp_x(n) = -disp_x(n) * thickness
        disp_y(n) = -disp_y(n) * thickness
        disp_z(n) = -disp_z(n) * thickness

        ! Far anchor nodes: prescribed zero displacement
        if(phi(n) .gt. 0.95) then  ! this 0.95 is a bit of ad-hoc
          fixed_nodes(n) = .true.
        end if

        ! Boundary-layer nodes: prescribed inward displacement
        if(mark_nodes(n) .gt. 0) then
          fixed_nodes(n) = .true.
        end if

      end do    ! through nodes

      !--------------------------!
      !   Estimate node normal   !
      !--------------------------!
      call Enlarge % Array_Real(node_nx, (/1, Grid % n_nodes/), 0.0)
      call Enlarge % Array_Real(node_ny, (/1, Grid % n_nodes/), 0.0)
      call Enlarge % Array_Real(node_nz, (/1, Grid % n_nodes/), 0.0)
      call Enlarge % Array_Int(node_cnt, (/1, Grid % n_nodes/), 0)

      do s = 1, Grid % n_faces
        if(Grid % faces_c(2,s) .lt. 0) then

          call Grid % Faces_Surface(s, sx, sy, sz)

          area = sqrt(sx*sx + sy*sy + sz*sz)
          Assert(area .gt. TINY)
          nx = sx / area
          ny = sy / area
          nz = sz / area

          do i_nod = 1, Grid % faces_n_nodes(s)
            n = Grid % faces_n(i_nod, s)
            node_nx(n)  = node_nx(n)  + nx
            node_ny(n)  = node_ny(n)  + ny
            node_nz(n)  = node_nz(n)  + nz
            node_cnt(n) = node_cnt(n) + 1
          end do
        end if
      end do

      do n = 1, Grid % n_nodes
        if(node_cnt(n) .gt. 0) then

          node_nx(n) = node_nx(n) / real(node_cnt(n))
          node_ny(n) = node_ny(n) / real(node_cnt(n))
          node_nz(n) = node_nz(n) / real(node_cnt(n))

          magn = sqrt(node_nx(n)**2 + node_ny(n)**2 + node_nz(n)**2)
          Assert(magn .gt. TINY)
          node_nx(n) = node_nx(n) / magn
          node_ny(n) = node_ny(n) / magn
          node_nz(n) = node_nz(n) / magn
        end if
      end do

      call Grid % Save_Debug_Vtu(append="node-normals",                        &
                                 vector_node = (/node_nx, node_ny, node_nz/),  &
                                 vector_name = "node-normals")

      !----------------------------------------!
      !   Call solver for node displacements   !
      !----------------------------------------!
      call Convert % Smooth_Node_Displacements(Grid,                       &
                                               fixed_nodes,                &
                                               disp_x,  disp_y,  disp_z,   &
                                               node_nx, node_ny, node_nz,  &
                                               1.0e-6, 1000)
      call Grid % Save_Vtu_Edges()

      deallocate(node_nx)
      deallocate(node_ny)
      deallocate(node_nz)
      deallocate(fixed_nodes)

      !--------------------------------------------------------!
      !   Leave the do i_lay loop after one compression pass   !
      !--------------------------------------------------------!
      exit

    end if  ! compress_domain

    !---------------------------------------------------------!
    !                                                         !
    !   If you are here, you are not compressing the domain   !
    !   but adding new layers (and nodes, faces and cells).   !
    !                                                         !
    !---------------------------------------------------------!

    ! Estimate to which new nodes will current nodes be projected
    cnt_n = 0
    call Enlarge % Array_Int(node_to, (/1, Grid % n_nodes/))
    do n = 1, Grid % n_nodes
      if(mark_nodes(n) .gt. 0) then
        cnt_n = cnt_n + 1
        node_to(n) = Grid % n_nodes + cnt_n
      end if
    end do

    !---------------------------------------------!
    !   Enlarge memory to accomodate new nodes    !
    !- - - - - - - - - - - - - - - - - - - - - - -!
    !   IMPORTANT: This changes Grid % n_nodes    !
    !---------------------------------------------!
    old_n = Grid % n_nodes
    call Grid % Allocate_Nodes(old_n + cnt_n)

    ! Set coordinates in these new nodes
    do n = 1, old_n
      if(mark_nodes(n) .gt. 0) then
        Grid % xn(node_to(n)) = Grid % xn(n) + disp_x(n) * thickness
        Grid % yn(node_to(n)) = Grid % yn(n) + disp_y(n) * thickness
        Grid % zn(node_to(n)) = Grid % zn(n) + disp_z(n) * thickness
      end if
    end do

    Assert(Grid % n_nodes .eq. old_n + cnt_n)

    deallocate(disp_x)
    deallocate(disp_y)
    deallocate(disp_z)

    !-------------------!
    !                   !
    !   Add new faces   !
    !                   !
    !-------------------!

    ! Estimate to which new faces will current faces be projected
    cnt_c = 0
    cnt_f = 0
    call Enlarge % Array_Int(face_to_face, (/1, Grid % n_faces/))
    call Enlarge % Array_Int(face_to_cell, (/1, Grid % n_faces/))
    do s = 1, Grid % n_faces
      if(mark_face(s) .eq. YES) then

        ! Handle cell mapping
        cnt_f = cnt_f + 1
        face_to_face(s) = Grid % n_faces + cnt_f

        ! Handle inside cell mapping too
        cnt_c = cnt_c + 1
        face_to_cell(s) = Grid % n_cells + cnt_c  ! new, additional cell
      end if
    end do
    Assert(cnt_f .eq. cnt_c)

    !----------------!
    !                !
    !   Find edges   !
    !                !
    !----------------!

    ! This are approximate sizes
    call Enlarge % Array_Int(edge_n1,  (/1, 4*cnt_f/))
    call Enlarge % Array_Int(edge_n2,  (/1, 4*cnt_f/))
    call Enlarge % Array_Int(edge_s1,  (/1, 4*cnt_f/))
    call Enlarge % Array_Int(edge_s2,  (/1, 4*cnt_f/))
    call Enlarge % Array_Int(edge_cnt, (/1, 4*cnt_f/))

    ! Counter for all edges, inner and outer
    cnt_e = 0

    !---------------------------------------------------------------!
    !                                                               !
    !   Find all edges touched by marked nodes on boundary faces.   !
    !   The two faces sharing each edge are stored in edge_s1 and   !
    !   edge_s2, so inner and outer side faces can be identified    !
    !   later from mark_face(s1) and mark_face(s2).                 !
    !                                                               !
    !---------------------------------------------------------------!
    do s = 1, Grid % n_faces

      c1 = Grid % faces_c(1, s)
      c2 = Grid % faces_c(2, s)

      !----------------------------------!
      !   Consider only boundary faces   !
      !----------------------------------!
      if(c2 .lt. 0) then

        ! Browse through face's nodes
        do i_nod = 1, Grid % faces_n_nodes(s)
          j_nod = i_nod + 1
          if(j_nod > Grid % faces_n_nodes(s)) j_nod = 1

          ni = Grid % faces_n(i_nod, s)
          nj = Grid % faces_n(j_nod, s)

          n1 = min(ni, nj)
          n2 = max(ni, nj)

          if(mark_nodes(n1) .gt. 0 .and. mark_nodes(n2) .gt. 0) then

            ! Check if this edge already exists ...
            e = 0
            do ee = 1, cnt_e
              if(edge_n1(ee) .eq. n1 .and. edge_n2(ee) .eq. n2) then
                e = ee
                exit
              end if
            end do

            ! ... and if not, add it to the list
            if(e == 0) then

              cnt_e = cnt_e + 1
              Assert(cnt_e .le. 4*cnt_f)

              ! Edge count will be needed to distinguish inner from outer edges
              edge_cnt(cnt_e) = edge_cnt(cnt_e) + 1

              edge_n1(cnt_e) = n1
              edge_n2(cnt_e) = n2

              ! Store first face touching this edge
              edge_s1(cnt_e) = s
              edge_s2(cnt_e) = 0

            ! Edge has been found
            else

              ! Edge count should be updated even if the edge was found
              ! (Actually, it is even more important in that case)
              edge_cnt(e) = edge_cnt(e) + 1

              ! Store the second face touching this edge
              Assert(edge_s1(e) .ne. 0)  ! should have been set already
              edge_s2(e) = s
            end if

          end if  ! two nodes (an edge) is marked
        end do    ! face's nodes
      end if      ! c2 .lt. 0
    end do        ! through faces

    ! All edges are counted twice now
    Assert(maxval(edge_cnt(1:cnt_e)) .eq. 2)
    Assert(minval(edge_cnt(1:cnt_e)) .eq. 2)

    ! De-allocate what you don't need
    deallocate(edge_cnt)

    !---------------------------------------------!
    !   Enlarge memory to accomodate new faces    !
    !- - - - - - - - - - - - - - - - - - - - - - -+-------------------!
    !   IMPORTANT: This changes Grid % n_faces and Grid % n_shadows   !
    !              although n_shadows is still zero at this point.    !
    !-----------------------------------------------------------------!
    old_f = Grid % n_faces
    Assert(Grid % n_shadows .eq. 0)  ! shadows shuld still be zero here

    call Allocate_Faces(Grid, Grid % n_faces + cnt_f + cnt_e, 0)
    Grid % n_faces = old_f

    ! Correct this: new grid will have some faces with four nodes for sure
    call Enlarge % Matrix_Int(Grid % faces_n,  &
                              i=(/1, 4/),      &
                              j=(/1,Grid % n_faces + cnt_f + cnt_e/))

    !---------------------------------------------------!
    !   Enlarge memory for cells (and boundary cells)   !
    !- - - - - - - - - - - - - - - - - - - - - - - - - -+---------------!
    !   IMPORTANT: This changes Grid % n_cells and Grid % n_bnd_cells   !
    !-------------------------------------------------------------------!
    old_c  = Grid % n_cells
    old_bc = Grid % n_bnd_cells
    call Grid % Allocate_Cells(Grid % n_cells     + cnt_c,   &
                               Grid % n_bnd_cells + cnt_e)
    Grid % n_cells     = old_c
    Grid % n_bnd_cells = old_bc

    !---------------------------------------------------------------!
    !                                                               !
    !   Create cells on top of formerly boundary faces              !
    !                                                               !
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
    !   For each marked boundary face, create one new cell using    !
    !   the face's old nodes and their projected counterparts.      !
    !   The old boundary face is then changed into an inside face   !
    !   by replacing its boundary cell with the new cell in         !
    !   faces_c(2,s).                                               !
    !                                                               !
    !   A new projected face is also created and connected to the   !
    !   new cell and to the old boundary cell, which keeps the      !
    !   original boundary condition at the projected surface.       !
    !---------------------------------------------------------------!
    do s = 1, Grid % n_faces
      if(mark_face(s) .eq. YES) then

        ! Copy nodes
        Grid % faces_n_nodes(face_to_face(s)) = Grid % faces_n_nodes(s)
        do i_nod = 1, Grid % faces_n_nodes(s)
          n = Grid % faces_n(i_nod, s)
          Grid % faces_n(i_nod, face_to_face(s)) = node_to(n)
        end do

        ! Copy cells
        c1 = Grid % faces_c(1, s)
        c2 = Grid % faces_c(2, s)
        Grid % faces_c(1, face_to_face(s)) = face_to_cell(s)
        Grid % faces_c(2, face_to_face(s)) = c2  ! keep the same boundary cell
        Grid % faces_c(2, s) = face_to_cell(s)

        Assert(face_to_cell(s) .gt. Grid % n_cells)

        !-------------------------------------------------!
        !   Triangular boundary face: create prism cell   !
        !-------------------------------------------------!
        if(Grid % faces_n_nodes(s) .eq. 3) then
          call Enlarge % Matrix_Int(Grid % cells_n, i=(/1,6/))
          Grid % cells_n_nodes(face_to_cell(s)) = 6  ! prism
          Grid % cells_n(1, face_to_cell(s)) = Grid % faces_n(1, s)
          Grid % cells_n(2, face_to_cell(s)) = Grid % faces_n(2, s)
          Grid % cells_n(3, face_to_cell(s)) = Grid % faces_n(3, s)
          Grid % cells_n(4, face_to_cell(s)) = node_to(Grid % faces_n(1, s))
          Grid % cells_n(5, face_to_cell(s)) = node_to(Grid % faces_n(2, s))
          Grid % cells_n(6, face_to_cell(s)) = node_to(Grid % faces_n(3, s))

        !---------------------------------------------------------!
        !   Quadrilateral boundary face: create hexahedral cell   !
        !---------------------------------------------------------!
        else if(Grid % faces_n_nodes(s) .eq. 4) then
          call Enlarge % Matrix_Int(Grid % cells_n, i=(/1,8/))
          Grid % cells_n_nodes(face_to_cell(s)) = 8  ! hexahedron
          Grid % cells_n(1, face_to_cell(s)) = Grid % faces_n(1, s)
          Grid % cells_n(2, face_to_cell(s)) = Grid % faces_n(2, s)
          Grid % cells_n(3, face_to_cell(s)) = Grid % faces_n(3, s)
          Grid % cells_n(4, face_to_cell(s)) = Grid % faces_n(4, s)
          Grid % cells_n(5, face_to_cell(s)) = node_to(Grid % faces_n(1, s))
          Grid % cells_n(6, face_to_cell(s)) = node_to(Grid % faces_n(2, s))
          Grid % cells_n(7, face_to_cell(s)) = node_to(Grid % faces_n(3, s))
          Grid % cells_n(8, face_to_cell(s)) = node_to(Grid % faces_n(4, s))

        !-----------------------------------------!
        !   Polygonal face: create a polyhedron   !
        !-----------------------------------------!
        else
          nn = Grid % faces_n_nodes(s)
          call Enlarge % Matrix_Int(Grid % cells_n, i=(/1,2*nn/))
          Grid % cells_n_nodes(face_to_cell(s)) = -2*nn  ! polyhedron
          do i_nod = 1, nn
            Grid % cells_n(i_nod,      face_to_cell(s)) = Grid % faces_n(i_nod, s)
            Grid % cells_n(i_nod + nn, face_to_cell(s)) = node_to(Grid % faces_n(i_nod, s))
          end do
        end if

        !-------------------!
        !   Boundary cell   !
        !-------------------!
        Assert(c2 .lt. 0)
        Assert(Grid%cells_n_nodes(c2) .eq. Grid%faces_n_nodes(face_to_face(s)))
        Grid % cells_n_nodes(c2) = Grid % faces_n_nodes(face_to_face(s))

        ! Make sure that boundary cells contains extruded
        ! nodes: copy them from Grid % faces_n structure
        do i_nod = 1, Grid % cells_n_nodes(c2)
          n = Grid % faces_n(i_nod, face_to_face(s))
          Assert(n .gt. 0)
          Grid % cells_n(i_nod, c2) = n
        end do

      end if  ! mark_face(s) .eq. YES
    end do    ! faces

    !-------------------------------------------------------------!
    !                                                             !
    !   Create new faces from the edges of the extruded regions   !
    !   For the outer edges, create new boundary cells as well.   !
    !                                                             !
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
    !   Do in two runs, first outer, then inner edges             !
    !     run = 1 -> skips inner edges, processes outer edges     !
    !     run = 2 -> skips outer edges, processes inner edges     !
    !-------------------------------------------------------------!
    cnt_ei = 0
    cnt_eb = 0
    do run = 1, 2

      do e = 1, cnt_e
        n1 = edge_n1(e)
        n2 = edge_n2(e)
        s1 = edge_s1(e)
        s2 = edge_s2(e)

        Assert(s1 .gt. 0)
        Assert(s2 .gt. 0)

        ! Is it inner our outer edge?
        inner =       mark_face(s1) .eq. YES  &
                .and. mark_face(s2) .eq. YES
        outer = .not. inner

        ! Run 1 creates side boundary faces, run 2 creates side inner faces.
        if(run .eq. 1 .and. inner) cycle
        if(run .eq. 2 .and. .not. inner) cycle

        ! At this moment, it is either:
        !   run==1 and .not. inner
        !   run==2 and inner

        if(outer) cnt_eb = cnt_eb + 1
        if(inner) cnt_ei = cnt_ei + 1

        !-------------------------------------------------------------!
        !   Set face nodes.  This is valid for both outer and inner   !
        !   faces.  Keep in mind that if you are in the outer loop    !
        !   (run == 1) the counter cnt_ei is 0                        !
        !-------------------------------------------------------------!

        ! Extruded face will be quadrilateral
        Grid % faces_n_nodes(Grid % n_faces + cnt_f + cnt_eb + cnt_ei) = 4

        ! Form the new face's nodes
        Grid % faces_n(1, Grid % n_faces + cnt_f + cnt_eb + cnt_ei) = n1
        Grid % faces_n(2, Grid % n_faces + cnt_f + cnt_eb + cnt_ei) = n2
        Grid % faces_n(3, Grid % n_faces + cnt_f + cnt_eb + cnt_ei) = node_to(n2)
        Grid % faces_n(4, Grid % n_faces + cnt_f + cnt_eb + cnt_ei) = node_to(n1)

        !---------------------------------!
        !   Set face cells in two runs:   !
        !     run == 1 outer              !
        !     run == 2 inner              !
        !---------------------------------!
        if(outer) then

          Assert(cnt_ei .eq. 0)

          ! Distinguish marked from boundary cell
          if(     mark_face(s1) .eq. YES  &
            .and. mark_face(s2) .eq. NO) then
            sm = s1  ! marked face
            sb = s2  ! non marked face
          else if(mark_face(s1) .eq. NO  &
            .and. mark_face(s2) .eq. YES) then
            sm = s2  ! marked face
            sb = s1  ! non-marked face
          else
            call Message % Error(60, "Something is wrong big time")
          end if

          ! Cell c1
          c1 = face_to_cell(sm)  ! new cell c1; projected from sm face
          Assert(c1 .gt. 0)
          Grid % faces_c(1, Grid % n_faces + cnt_f + cnt_eb) = c1

          ! Cell c2
          Assert(Grid % faces_c(2, sb) .lt. 0)
          c2 = -Grid % n_bnd_cells - cnt_eb  ! new cell c2

          ! Copy the old region information
          Grid % region % at_cell(c2)  &
            = Grid % region % at_cell(Grid % faces_c(2, sb))
          Assert(Grid % region % at_cell(c2) .gt. 0)

          ! Set c2
          Grid % faces_c(2, Grid % n_faces + cnt_f + cnt_eb) = c2

          ! Copy nodes from the newly formed face to the new boundary cell
          Grid % cells_n_nodes(c2)  &
            = Grid % faces_n_nodes(Grid % n_faces + cnt_f + cnt_eb)
          Assert(Grid % cells_n_nodes(c2) .eq. 4)
          Grid % cells_n(1:4, c2)  &
            = Grid % faces_n(1:4, Grid % n_faces + cnt_f + cnt_eb)

        !-------------------------------------------!
        !   In the second run, do the inner edges   !
        !-------------------------------------------!
        else if(inner) then

          c1 = min(face_to_cell(s1), face_to_cell(s2))
          c2 = max(face_to_cell(s1), face_to_cell(s2))

          Assert(c1 .gt. Grid % n_cells)
          Assert(c2 .gt. Grid % n_cells)

          Grid % faces_c(1, Grid % n_faces + cnt_f + cnt_eb + cnt_ei) = c1
          Grid % faces_c(2, Grid % n_faces + cnt_f + cnt_eb + cnt_ei) = c2
        end if

      end do

    end do      ! run
    Assert(cnt_e .eq. cnt_ei + cnt_eb)

    ! Print more info (spaces here are aligned with print statements above)
    print *, "# Additional boundary faces:       ", cnt_eb
    print *, "# Additional inside faces:         ", cnt_ei
    print *, "# Total additional faces:          ", cnt_e

    deallocate(edge_n1)
    deallocate(edge_n2)
    deallocate(edge_s1)
    deallocate(edge_s2)

    deallocate(mark_nodes)
    deallocate(mark_face)
    deallocate(node_to)
    deallocate(face_to_face)
    deallocate(face_to_cell)

    !-----------------------------------------------------------------------!
    !                                                                       !
    !   IMPORTANT: Increase the number of faces, cells and boundary cells   !
    !                                                                       !
    !-----------------------------------------------------------------------!
    Grid % n_faces     = Grid % n_faces     + cnt_f + cnt_e
    Grid % n_cells     = Grid % n_cells     + cnt_c
    Grid % n_bnd_cells = Grid % n_bnd_cells + cnt_eb

    !-----------------------------------------------------------------------!
    !                                                                       !
    !   Sort new faces so that they resemble the order in which they were   !
    !   before inserting the boundary layer.  That means, boundary faces    !
    !   first, followed by inside faces, each sorted by c1 index.           !
    !                                                                       !
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
    !   This is probably not very important, but keeping a familiar order   !
    !   makes debugging and comparison of converted grids much easier.      !
    !-----------------------------------------------------------------------!
    call Enlarge % Array_Int(key, i=(/1, Grid % n_faces/))
    do s = 1, Grid % n_faces
      c1 = Grid % faces_c(1, s)
      c2 = Grid % faces_c(2, s)

      ! The key is formed so that the faces are ordered as:
      !   1. boundary faces first
      !   2. inside faces after boundary faces
      !   3. within each group, increasing c1
      !   4. for equal c1, preserve old face order via + s
      if(c2 .lt. 0) then
        key(s) = c1 * (Grid % n_faces + 1) + s
      else
        key(s) = (Grid % n_cells + c1) * (Grid % n_faces + 1) + s
      end if

      Grid % old_f(s) = s
    end do

    call Sort % Int_Carry_Int(key, Grid % old_f)

    do s = 1, Grid % n_faces
      Grid % new_f(Grid % old_f(s)) = s
    end do

    call Grid % Sort_Faces_By_Index(Grid % new_f)

    !----------------------------!
    !                            !
    !   Repair face directions   !
    !                            !
    !----------------------------!
    do run = 1, 2  ! correct, then check
      do s = 1, Grid % n_faces

        c1 = Grid % faces_c(1, s)
        c2 = Grid % faces_c(2, s)
        call Grid % Faces_Surface(s, sx, sy, sz)
        call Grid % Cells_Center(c1, xc1, yc1, zc1)
        call Grid % Cells_Center(c2, xc2, yc2, zc2)
        dx = xc2 - xc1
        dy = yc2 - yc1
        dz = zc2 - zc1
        dot = sx*dx + sy*dy + sz * dz

        !------------------------------------------------------------!
        !   In the first run, the order of nodes is not guaranteed   !
        !------------------------------------------------------------!
        if(run .eq. 1 .and. dot .lt. 0.0) then

          ! Triangle: [1,2,3] -> [1,3,2]
          if(Grid % faces_n_nodes(s) .eq. 3) then
            call Swap_Int(Grid % faces_n(2, s), Grid % faces_n(3, s))

          ! Quad: [1,2,3,4] -> [1,4,3,2]
          else if(Grid % faces_n_nodes(s) .eq. 4) then
            call Swap_Int(Grid % faces_n(2, s), Grid % faces_n(4, s))

          ! Polygons: [1,2,3,4,5]     -> [1,5,4,3,2]
          !           [1,2,3,4,5,6]   -> [1,6,5,4,3,2]
          !           [1,2,3,4,5,6,7] -> [1,7,6,5,4,3,2]
          else
            call Sort % Reverse_Order_Int(  &
              Grid % faces_n(2:Grid % faces_n_nodes(s), s))
          end if
        end if

        !--------------------------------------------------!
        !   In the second run, all faces should be fixed   !
        !--------------------------------------------------!
        if(run .eq. 2) then
          Assert(dot .ge. 0)
        end if

      end do  ! faces
    end do    ! run

    !------------------------------!
    !                              !
    !   Check as much as you can   !
    !                              !
    !------------------------------!
    do s = 1, Grid % n_faces
      c1 = Grid % faces_c(1, s)
      c2 = Grid % faces_c(2, s)
      Assert(c1 .gt. 0)
      Assert(c2 .ne. 0)
      Assert(c1 .ne. c2)
      if(c2 .gt. 0) then
        Assert(c1 .lt. c2)
      end if

      Assert(Grid % faces_n_nodes(s) .gt. 0)
      do i_nod = 1, Grid % faces_n_nodes(s)
        n = Grid % faces_n(i_nod, s)
        Assert(n .gt. 0)
        nn1 = abs(Grid % cells_n_nodes(c1))
        nn2 = abs(Grid % cells_n_nodes(c2))
        Assert(any(Grid % cells_n(1:nn1, c1) .eq. n))
        Assert(any(Grid % cells_n(1:nn2, c2) .eq. n))
      end do
    end do

    !------------------------!
    !                        !
    !   Save for debugging   !
    !                        !
    !------------------------!
    if(DEBUG) then
    !@  Grid % s(:) = 1.0
    !@  call Grid % Save_Vtu_Faces(sub=(/0,0/), volume_flux=mark_face)
    !@  call Enlarge % Array_Real(mark_nodes, (/1, Grid % n_nodes/))
    !@
    !@  call Grid % Save_Debug_Vtu(append="node-count",         &
    !@                             scalar_node = mark_nodes,    &
    !@                             scalar_name = "node-count")
    !@  call Enlarge % Array_Real(disp_x, (/1, Grid % n_nodes/))
    !@  call Enlarge % Array_Real(disp_y, (/1, Grid % n_nodes/))
    !@  call Enlarge % Array_Real(disp_z, (/1, Grid % n_nodes/))
    !@  call Grid % Save_Debug_Vtu(append="node-thicknesss",    &
    !@                             vector_node = (/disp_x,disp_y,disp_z/),  &
    !@                             vector_name = "node-thickness")
    end if

  end do  ! through layers

  end subroutine
