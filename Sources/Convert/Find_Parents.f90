!==============================================================================!
  subroutine Find_Parents(grid)
!------------------------------------------------------------------------------!
!   Looks boundary cells' parents for meshes in which they are not given       !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer              :: fn(6,4), i_fac, j_fac, i_nod, i_cel, c1, c2
  integer              :: nodes(4), n_match, dir, bc, cnt_c, cnt_f
  integer              :: n_face_nodes ! number of nodes in a face
  integer              :: n_cell_faces ! number of faces in a cell
  logical, allocatable :: is_node_bnd(:)
  integer, allocatable :: cell_near_bnd(:)
  integer, allocatable :: cr1(:), cr2(:), cr3(:), w1(:), w2(:)  ! for sorting
!------------------------------------------------------------------------------!
  include 'Cells_Faces_Nodes.f90'
!==============================================================================!

  print *, '#================================================='
  print *, '# Parent information not given in the input file!'
  print *, '# Looking for parents. This may take a few minutes'
  print *, '#-------------------------------------------------'

  print *, '# Number of boundary cells: ', grid % n_bnd_cells
  print *, '# Number of inside cells:   ', grid % n_cells

  if(grid % n_cells .eq. 0) then
    print *, '# Number of cells inside the domain is zero'
    print *, '# Are you sure you meshed the domain in 3D?'
    print *, '# Exiting!'
    print *, '#-------------------------------------------------'
    stop
  end if

  allocate(is_node_bnd(grid % n_nodes));    is_node_bnd(:)   = .false.
  allocate(cell_near_bnd(grid % n_cells));  cell_near_bnd(:) = 0

  !--------------------------------------------------------------!
  !   Find all cells near boundary and allocate working arrays   !
  !--------------------------------------------------------------!

  ! Mark all nodes on the boundary
  is_node_bnd(:) = .false.
  do c2 = -grid % n_bnd_cells, -1
    if( grid % bnd_cond % color(c2) .gt. 0 ) then
      do i_nod = 1, grid % cells_n_nodes(c2)  ! 3 or 4
        is_node_bnd( grid % cells_n(i_nod, c2) ) = .true.
      end do
    end if
  end do

  ! Count cells near boundaries
  cnt_c = 0
  cnt_f = grid % n_bnd_cells
  do c1 = 1, grid % n_cells
    do i_nod = 1, grid % cells_n_nodes(c1)  ! 4 to 8
      if( is_node_bnd(grid % cells_n(i_nod, c1)) ) then
        cnt_c = cnt_c + 1
        cell_near_bnd(cnt_c) = c1
        if(grid % cells_n_nodes(c1) .eq. 4) cnt_f = cnt_f + 4  ! tet
        if(grid % cells_n_nodes(c1) .eq. 5) cnt_f = cnt_f + 5  ! pyr
        if(grid % cells_n_nodes(c1) .eq. 6) cnt_f = cnt_f + 5  ! wed
        if(grid % cells_n_nodes(c1) .eq. 8) cnt_f = cnt_f + 6  ! hex
        goto 1
      end if
    end do
1   continue
  end do
  PRINT *, '# Number of all cells near boundary:          ', cnt_c
  PRINT *, '# Number of all cells''s faces near boundary: ', cnt_f

  allocate(cr1(cnt_f));  cr1(:) = 0
  allocate(cr2(cnt_f));  cr2(:) = 0
  allocate(cr3(cnt_f));  cr3(:) = 0
  allocate(w1 (cnt_f));  w1 (:) = 0
  allocate(w2 (cnt_f));  w2 (:) = 0

  !----------------------!
  !   Real work begins   !
  !----------------------!

  do bc = 1, grid % n_bnd_cond

    cr1(:) = HUGE_INT
    cr2(:) = HUGE_INT
    cr3(:) = HUGE_INT
    w1 (:) = HUGE_INT
    w2 (:) = HUGE_INT

    cnt_f = 0
    do c2 = -grid % n_bnd_cells, -1
      if( grid % bnd_cond % color(c2) .eq. bc ) then

        ! Increase total face count
        cnt_f = cnt_f + 1

        ! Form and store local nodes list
        nodes(:) = HUGE_INT
        do i_nod = 1, grid % cells_n_nodes(c2)  ! 3 or 4
          nodes(i_nod) = grid % cells_n(i_nod, c2)
        end do
        call Sort % Int_Array(nodes(1:4))
        cr1(cnt_f) = nodes(1)
        cr2(cnt_f) = nodes(2)
        cr3(cnt_f) = nodes(3)
        w1 (cnt_f) = c2        ! not sure this will be needed, but OK
        w2 (cnt_f) = bc        ! not sure this will be needed, but OK
      end if
    end do

    do i_cel = 1, cnt_c
      c1 = cell_near_bnd(i_cel)

      ! Fetch the faces from this cell
      if(grid % cells_n_nodes(c1) .eq. 4) fn = tet
      if(grid % cells_n_nodes(c1) .eq. 5) fn = pyr
      if(grid % cells_n_nodes(c1) .eq. 6) fn = wed
      if(grid % cells_n_nodes(c1) .eq. 8) fn = hex

      n_cell_faces = 6  ! assume it is a hexahedron
      if(grid % cells_n_nodes(c1) .eq. 4) n_cell_faces = 4
      if(grid % cells_n_nodes(c1) .eq. 5) n_cell_faces = 5
      if(grid % cells_n_nodes(c1) .eq. 6) n_cell_faces = 5

      ! Browse through all possible faces and their nodes
      do i_fac = 1, n_cell_faces

        ! Increase total face count
        cnt_f = cnt_f + 1

        ! Estimate number of nodes for this face
        n_face_nodes = 4             ! assume it is a quad
        if( fn(i_fac, 4) < 0 ) then  ! nope, it was a triangle
          n_face_nodes = 3
        end if

        ! Form and store local nodes list
        nodes(:) = HUGE_INT
        do i_nod = 1, n_face_nodes
          if(fn(i_fac, i_nod) > 0) then  ! if this node exists
            nodes(i_nod) = grid % cells_n(fn(i_fac, i_nod), c1)
          end if
        end do
        call Sort % Int_Array(nodes(1:4))
        cr1(cnt_f) = nodes(1)
        cr2(cnt_f) = nodes(2)
        cr3(cnt_f) = nodes(3)
        w1 (cnt_f) = c1
        w2 (cnt_f) = i_fac

      end do  ! i_fac

    end do  ! c1

    call Sort % Three_Int_Carry_Two_Int(                   &
                cr1(1:cnt_f), cr2(1:cnt_f), cr3(1:cnt_f),  &
                w1(1:cnt_f), w2(1:cnt_f))

    n_match = 0
    do i_fac = 1, cnt_f-1
      j_fac = i_fac + 1
      if( w1(i_fac) * w1(j_fac) < 0 ) then  ! they have different signs

        ! Is it a match?
        if( cr1(i_fac) .eq. cr1(j_fac) .and.  &
            cr2(i_fac) .eq. cr2(j_fac) .and.  &
            cr3(i_fac) .eq. cr3(j_fac) ) then

          n_match = n_match + 1
          if( w1(i_fac) > 0 ) then  ! i_fac stores internal cell
            c1  = w1(i_fac)
            dir = w2(i_fac)
          else                      ! j_fac stores internal cell
            c1  = w1(j_fac)
            dir = w2(j_fac)
          end if
          grid % cells_bnd_color(dir, c1) = bc

        end if
      end if
    end do
    print *, '# Number of matches: ', n_match

  end do  ! bc

  end subroutine
