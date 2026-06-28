!==============================================================================!
  subroutine Insert_Layers(Convert, Grid)
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
  class(Convert_Type)        :: Convert  !! parent class
  type(Grid_Type),    target :: Grid     !! primal grid
!------------------------------[Local parameters]------------------------------!
  logical, parameter :: DEBUG = .false.
!-----------------------------------[Locals]-----------------------------------!
  character(SL)        :: answer
  integer              :: n_inserting_regions, i_reg, reg, s, c1, c2, ni, nj
  integer              :: cnt_c, cnt_f, cnt_n, cnt_e, cnt_ei, cnt_eb, run
  integer              :: n, i_nod, j_nod, e, ee, n1, n2, s1, s2, sm, sb
  integer              :: old_n, old_f, old_c, old_bc
  real                 :: area, dot, eps
  real                 :: fnx, fny, fnz
  real                 :: xf, yf, zf, xc, yc, zc, sx, sy, sz
  real                 :: xc1, yc1, zc1, xc2, yc2, zc2, dx, dy, dz
  real                 :: shift, a(3,3), b(3), x(3)
  logical              :: invertible, inner, outer, warning
  logical, allocatable :: edge_on_bnd(:)  ! is edge on the boundary
  integer, allocatable :: node_to(:), face_to_face(:), face_to_cell(:)
  integer, allocatable :: nodes_rank(:), edge_cnt(:), key(:)
  integer, allocatable :: edge_n1(:), edge_n2(:), edge_s1(:), edge_s2(:)
  real,    allocatable :: mark_nodes(:), mark_faces(:)
  real,    allocatable :: nx(:), ny(:), nz(:)
  real,    allocatable :: a11(:), a12(:), a13(:)
  real,    allocatable :: a22(:), a23(:), a33(:)
  real,    allocatable :: bx(:),  by(:),  bz(:)
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Convert)
!==============================================================================!

  do

    call Print_Regions_List(Grid)

    call Message % Framed(60,                                              &
      "Inserting boundary layers.                                    ",    &
      "Type ordinal number(s) of boundary-condition regions to       " //  &
      "expand, followed by the boundary-layer thickness.             " //  &
      "\n \n                                                         " //  &
      "Example:  2  5  7  0.01                                       " //  &
      "expands regions 2, 5 and 7 by a layer of thickness 0.01.      " //  &
      "\n \n                                                         " //  &
      "Type skip to skip insertion of boundary layers")

    call File % Read_Line(5)
    answer = Line % tokens(1)
    call String % To_Upper_Case(answer)

    if(answer .eq. 'SKIP') then
      print *, "# Finished inserting boundary layers."
      exit
    end if

    !-----------------------------------------------!
    !                                               !
    !   Users wants to introduce a boundary layer   !
    !                                               !
    !-----------------------------------------------!

    ! Just to see how is the grid sorted before
    !@ do s = 1, Grid % n_faces
    !@   c1 = Grid % faces_c(1,s)
    !@   c2 = Grid % faces_c(2,s)
    !@   if(c2 .gt. 0) then
    !@     write(8, *) s, c1, c2
    !@   else
    !@     write(8, *) s, c1, c2, Grid % region % at_cell(c2)
    !@   end if
    !@ end do

    print *, '# Inserting boundary layer!'

    !------------------------------------------------!
    !   Fetch the number of regions to be extruded   !
    !------------------------------------------------!
    n_inserting_regions = Line % n_tokens - 1  ! the last one is shift

    !---------------------!
    !   Fetch the shift   !
    !---------------------!
    read(Line % tokens(Line % n_tokens), *)  shift
    print '(a,es10.3)', " # Inserting a boundary layer with thicness ", shift

    !------------------------------------------------------!
    !                                                      !
    !   Sort nodes by boundary regions you are extruding   !
    !       in order to save memory on node matrices       !
    !                                                      !
    !------------------------------------------------------!
    call Enlarge % Array_Int(nodes_rank, (/1, Grid % n_nodes/))
    nodes_rank(:) = HUGE_INT

    cnt_n = 0  ! number of extruded nodes
    do i_reg = 1, n_inserting_regions
      read(Line % tokens(i_reg), *) reg

      ! Browse through all faces
      do s = 1, Grid % n_faces
        c2 = Grid % faces_c(2,s)

        ! Check if face belings to the selected region ...
        if(c2 < 0) then
          if(Grid % region % at_cell(c2) .eq. reg) then

            ! ... and if so, mark the face's nodes
            do i_nod = 1, Grid % faces_n_nodes(s)
              n = Grid % faces_n(i_nod, s)
              if(nodes_rank(n) .eq. HUGE_INT) then
                cnt_n = cnt_n + 1
                nodes_rank(n) = cnt_n  ! extruded node rank
              end if
            end do
          end if
        end if
      end do  ! through faces
    end do    ! through regions

    ! The sorting can take place now
    call Grid % Sort_Nodes_By_Index(nodes_rank)

    print *, "# Old number of nodes:             ", Grid % n_nodes
    print *, "# New nodes in the boundary layer: ", cnt_n

    call Enlarge % Array_Real(mark_nodes, (/1, Grid % n_nodes/))
    call Enlarge % Array_Real(mark_faces, (/1, Grid % n_faces/))
    call Enlarge % Array_Real(nx,         (/1, Grid % n_nodes/))
    call Enlarge % Array_Real(ny,         (/1, Grid % n_nodes/))
    call Enlarge % Array_Real(nz,         (/1, Grid % n_nodes/))

    call Enlarge % Array_Real(a11, (/1, cnt_n/))
    call Enlarge % Array_Real(a12, (/1, cnt_n/))
    call Enlarge % Array_Real(a13, (/1, cnt_n/))
    call Enlarge % Array_Real(a22, (/1, cnt_n/))
    call Enlarge % Array_Real(a23, (/1, cnt_n/))
    call Enlarge % Array_Real(a33, (/1, cnt_n/))
    call Enlarge % Array_Real(bx,  (/1, cnt_n/))
    call Enlarge % Array_Real(by,  (/1, cnt_n/))
    call Enlarge % Array_Real(bz,  (/1, cnt_n/))

    mark_nodes(:) = 0
    mark_faces(:) = 0

    nx (:) = 0.0;  ny (:) = 0.0;  nz (:) = 0.0
    a11(:) = 0.0;  a12(:) = 0.0;  a13(:) = 0.0
    a22(:) = 0.0;  a23(:) = 0.0;  a33(:) = 0.0
    bx (:) = 0.0;  by (:) = 0.0;  bz (:) = 0.0

    warning = .false.
    do i_reg = 1, n_inserting_regions
      read(Line % tokens(i_reg), *) reg

      do s = 1, Grid % n_faces
        c1 = Grid % faces_c(1,s)
        c2 = Grid % faces_c(2,s)

        if(c2 < 0) then
          if(Grid % region % at_cell(c2) .eq. reg) then

            if(mark_faces(s) .eq. 0) then
              mark_faces(s) = 1

              ! Layer can be added to primal grids only, meaning the
              ! number of faces's nodes can be only three or four
              Assert(Grid % faces_n_nodes(s) .ge. 3)
              Assert(Grid % faces_n_nodes(s) .le. 4)

              ! Calculate face's surface normal
              call Grid % Faces_Surface(s, sx, sy, sz)

              area = sqrt(sx**2 + sy**2 + sz**2)
              Assert(area > TINY)

              fnx = sx / area
              fny = sy / area
              fnz = sz / area

              call Grid % Faces_Center(s, xf, yf, zf)

              call Grid % Cells_Center(c1, xc, yc, zc)

              dot = fnx*(xf-xc) + fny*(yf-yc) + fnz*(zf-zc)
              if(dot < 0.0) then
                warning = .true.
                fnx = -fnx
                fny = -fny
                fnz = -fnz
              end if

              do i_nod = 1, Grid % faces_n_nodes(s)
                n = Grid % faces_n(i_nod, s)
                Assert(n .gt. 0)
                Assert(n .le. cnt_n)

                mark_nodes(n) = mark_nodes(n) + 1

                a11(nodes_rank(n)) = a11(nodes_rank(n)) + fnx * fnx
                a12(nodes_rank(n)) = a12(nodes_rank(n)) + fnx * fny
                a13(nodes_rank(n)) = a13(nodes_rank(n)) + fnx * fnz
                a22(nodes_rank(n)) = a22(nodes_rank(n)) + fny * fny
                a23(nodes_rank(n)) = a23(nodes_rank(n)) + fny * fnz
                a33(nodes_rank(n)) = a33(nodes_rank(n)) + fnz * fnz

                bx(nodes_rank(n)) = bx(nodes_rank(n)) + fnx
                by(nodes_rank(n)) = by(nodes_rank(n)) + fny
                bz(nodes_rank(n)) = bz(nodes_rank(n)) + fnz
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

    eps = MICRO

    do n = 1, Grid % n_nodes
      if(mark_nodes(n) .gt. 0) then

        a(1,1) = a11(nodes_rank(n)) + eps
        a(1,2) = a12(nodes_rank(n))
        a(1,3) = a13(nodes_rank(n))
        a(2,1) = a12(nodes_rank(n))
        a(2,2) = a22(nodes_rank(n)) + eps
        a(2,3) = a23(nodes_rank(n))
        a(3,1) = a13(nodes_rank(n))
        a(3,2) = a23(nodes_rank(n))
        a(3,3) = a33(nodes_rank(n)) + eps

        b(1) = bx(nodes_rank(n))
        b(2) = by(nodes_rank(n))
        b(3) = bz(nodes_rank(n))
        x(:) = 0.0

        call Math % Gaussian_Elimination(3, a, b, x, invertible)

        if(invertible) then
          nx(n) = x(1)
          ny(n) = x(2)
          nz(n) = x(3)
        else
          nx(n) = bx(nodes_rank(n)) / mark_nodes(n)
          ny(n) = by(nodes_rank(n)) / mark_nodes(n)
          nz(n) = bz(nodes_rank(n)) / mark_nodes(n)
        end if

      end if  ! mark_nodes(n) .gt. 0
    end do    ! through nodes

    deallocate(nodes_rank)
    deallocate(a11)
    deallocate(a12)
    deallocate(a13)
    deallocate(a22)
    deallocate(a23)
    deallocate(a33)

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
      if(mark_nodes(n) > 0.0) then
        Grid % xn(node_to(n)) = Grid % xn(n) + nx(n) * shift
        Grid % yn(node_to(n)) = Grid % yn(n) + ny(n) * shift
        Grid % zn(node_to(n)) = Grid % zn(n) + nz(n) * shift
      end if
    end do

    Assert(Grid % n_nodes .eq. old_n + cnt_n)

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
      if(mark_faces(s) .gt. 0) then

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
    call Enlarge % Array_Int(edge_n1,     (/1, 3*cnt_f/))
    call Enlarge % Array_Int(edge_n2,     (/1, 3*cnt_f/))
    call Enlarge % Array_Int(edge_s1,     (/1, 3*cnt_f/))
    call Enlarge % Array_Int(edge_s2,     (/1, 3*cnt_f/))
    call Enlarge % Array_Int(edge_cnt,    (/1, 3*cnt_f/))
    call Enlarge % Array_Log(edge_on_bnd, (/1, 3*cnt_f/))

    cnt_e = 0

    !---------------------------------------------------------------------!
    !                                                                     !
    !   This will find inner and outer edges, distinguished by edge_cnt   !
    !                                                                     !
    !---------------------------------------------------------------------!
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
              Assert(cnt_e .le. 3*cnt_f)

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
      if(mark_faces(s) > 0.0) then

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

        if(Grid % faces_n_nodes(s) .eq. 3) then
          call Enlarge % Matrix_Int(Grid % cells_n, i=(/1,6/))
          Grid % cells_n_nodes(face_to_cell(s)) = 6  ! prism
          Grid % cells_n(1, face_to_cell(s)) = Grid % faces_n(1, s)
          Grid % cells_n(2, face_to_cell(s)) = Grid % faces_n(2, s)
          Grid % cells_n(3, face_to_cell(s)) = Grid % faces_n(3, s)
          Grid % cells_n(4, face_to_cell(s)) = node_to(Grid % faces_n(1, s))
          Grid % cells_n(5, face_to_cell(s)) = node_to(Grid % faces_n(2, s))
          Grid % cells_n(6, face_to_cell(s)) = node_to(Grid % faces_n(3, s))

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
        end if

        ! Take care of the boundary cell
        Assert(c2 .lt. 0)

        Assert(Grid % cells_n_nodes(c2) .eq. Grid % faces_n_nodes(face_to_face(s)))
        Grid % cells_n_nodes(c2) = Grid % faces_n_nodes(face_to_face(s))

        do i_nod = 1, Grid % cells_n_nodes(c2)
          n = Grid % faces_n(i_nod, face_to_face(s))
          Assert(n .gt. 0)
          Grid % cells_n(i_nod, c2) = n
        end do

      end if  ! mark_faces(s) .gt. 0
    end do    ! faces

    !---------------------------------------------------------!
    !   Do in two runs, first outer, then inner edges         !
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
    !   run = 1 -> skips inner edges, processes outer edges   !
    !   run = 2 -> skips outer edges, processes inner edges   !
    !---------------------------------------------------------!
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
        inner = mark_faces(s1) .gt. 0 .and. mark_faces(s2) .gt. 0
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
          if(     mark_faces(s1) .gt. 0  &
            .and. mark_faces(s2) .eq. 0) then
            sm = s1  ! marked face
            sb = s2  ! non marked face
          else if(mark_faces(s1) .eq. 0  &
            .and. mark_faces(s2) .gt. 0) then
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

    !---------------------------------------------!
    !   IMPORTANT: Increase the number of faces   !
    !---------------------------------------------!
    Grid % n_faces = Grid % n_faces + cnt_f + cnt_e

    !---------------------------------------------!
    !   IMPORTANT: Increase the number of cells   !
    !---------------------------------------------!
    Grid % n_cells     = Grid % n_cells     + cnt_c
    Grid % n_bnd_cells = Grid % n_bnd_cells + cnt_eb

    !-----------------------------------------------------------------------!
    !                                                                       !
    !   Sort new faces so that they resemble the order in which they were   !
    !   before inserting the boundary layer.  That means, boundary faces    !
    !   first, followed by inside faces, each sorted by c1 index.           !
    !                                                                       !
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
    !   This is probably not very important, but I was getting stressed     !
    !   with the new oder, or lacko f, of the faces in the new grid.        !
    !-----------------------------------------------------------------------!
    call Enlarge % Array_Int(key, i=(/1, Grid % n_faces/))
    do s = 1, Grid % n_faces
      c1 = Grid % faces_c(1, s)
      c2 = Grid % faces_c(2, s)

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

    !--------------------------------!
    !                                !
    !   Repair new face directions   !
    !                                !
    !--------------------------------!
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
      if(run .eq. 1 .and. dot .lt. 0.0) then
        if(Grid % faces_n_nodes(s) .eq. 3) then
          call Swap_Int(Grid % faces_n(2, s), Grid % faces_n(3, s))
        else if(Grid % faces_n_nodes(s) .eq. 4) then
          call Swap_Int(Grid % faces_n(2, s), Grid % faces_n(4, s))
        else
          call Message % Error(60,                                      &
            "Something is wrong big time: a face which has neither "//  &
            "three nor four nodes have been detected",                  &
            file= __FILE__, line = __LINE__)
        end if
      end if
      if(run .eq. 2) then
        Assert(dot .ge. 0)
      end if
    end do
    end do

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
        Assert(any(Grid % cells_n(1:Grid % cells_n_nodes(c1), c1) .eq. n))
        Assert(any(Grid % cells_n(1:Grid % cells_n_nodes(c2), c2) .eq. n))
      end do
    end do

    ! Just to see how is the grid sorted after
    !@ do s = 1, Grid % n_faces
    !@   c1 = Grid % faces_c(1,s)
    !@   c2 = Grid % faces_c(2,s)
    !@   if(c2 .gt. 0) then
    !@     write(9, *) s, c1, c2
    !@   else
    !@     write(9, *) s, c1, c2, Grid % region % at_cell(c2)
    !@   end if
    !@ end do

    !------------------------!
    !                        !
    !   Save for debugging   !
    !                        !
    !------------------------!
      Grid % s(:) = 1.0
      call Grid % Save_Vtu_Faces(sub=(/0,0/), volume_flux=mark_faces)
    if(DEBUG) then
      Grid % s(:) = 1.0
      call Grid % Save_Vtu_Faces(sub=(/0,0/), volume_flux=mark_faces)
      call Enlarge % Array_Real(mark_nodes, (/1, Grid % n_nodes/))

      call Grid % Save_Debug_Vtu(append="node-count",         &
                                 scalar_node = mark_nodes,    &
                                 scalar_name = "node-count")
      call Enlarge % Array_Real(nx,         (/1, Grid % n_nodes/))
      call Enlarge % Array_Real(ny,         (/1, Grid % n_nodes/))
      call Enlarge % Array_Real(nz,         (/1, Grid % n_nodes/))
      call Grid % Save_Debug_Vtu(append="node-shifts",        &
                                 vector_node = (/nx,ny,nz/),  &
                                 vector_name = "node-shifts")
      STOP
    end if
  end do

  end subroutine
