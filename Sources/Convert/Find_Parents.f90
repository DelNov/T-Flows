!==============================================================================!
  subroutine Find_Parents(grid)
!------------------------------------------------------------------------------!
!   Looks boundary cells' parents for meshes in which they are not given       !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod, only: HUGE_INT
  use Grid_Mod,  only: Grid_Type
  use Sort_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!------------------------------------------------------------------------------!
  include 'Cell_Numbering_Neu.f90'
!-----------------------------------[Locals]-----------------------------------!
  integer              :: fn(6,4), j, n1, n2, c1, c2, cb, run
  integer              :: n_match
  integer              :: n_bnd_nodes  ! near boundary nodes
  integer              :: n_bnd_proj   ! bnd. cells projected from cells
  integer              :: n_near_bnd   ! number of near boundary cells
  integer              :: n_face_nodes ! number of nodes in a face
  integer              :: n_cell_faces ! number of faces in a cell
  integer              :: n_fr_cells   ! for statistics
  integer, allocatable :: near_bnd(:)  ! near boundary cells
  logical, allocatable :: is_node_bnd(:)
  logical, allocatable :: is_cell_marked(:)
  integer, allocatable :: a1(:), a2(:), a3(:), b(:), c(:), d(:)
  integer              :: ai(4)
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

  !----------------------------------------!
  !   Allocate memory for helping arrays   !
  !----------------------------------------!

  ! Needed for sorting
  allocate(a1(grid % n_bnd_cells))
  allocate(a2(grid % n_bnd_cells))
  allocate(a3(grid % n_bnd_cells))
  allocate(b (grid % n_bnd_cells))
  allocate(c (grid % n_bnd_cells))
  allocate(d (grid % n_bnd_cells))

  ! Needed for algorithm logic
  allocate(is_node_bnd   (grid % n_nodes))
  allocate(is_cell_marked(grid % n_cells))
  is_node_bnd   (:) = .false.
  is_cell_marked(:) = .false.

  !----------------------------------------------------------!
  !   Prepare near boundary cells to work on a smaller set   !
  !----------------------------------------------------------!

  ! Mark all boundary nodes
  do c2 = -grid % n_bnd_cells, -1
    do n2 = 1, grid % cells_n_nodes(c2)  ! 3 or 4
      is_node_bnd( grid % cells_n(n2, c2) ) = .true.
    end do
  end do

  ! Count boundary nodes (just for information)
  n_bnd_nodes = 0
  do n1 = 1, grid % n_nodes
    if( is_node_bnd(n1) ) n_bnd_nodes = n_bnd_nodes + 1
  end do
  print *, '# Total    nodes: ', grid % n_nodes
  print *, '# Boundary nodes: ', n_bnd_nodes

  ! Store near boundary cells
  do run = 1, 2
    n_near_bnd        = 0
    is_cell_marked(:) = .false.
    do c1 = 1, grid % n_cells

      do n1 = 1, grid % cells_n_nodes(c1)  ! 4 to 8

        ! If at least one node is on the boundary, cell is there too
        if( is_node_bnd( grid % cells_n(n1, c1) ) ) then
          n_near_bnd = n_near_bnd + 1
          if(run .eq. 2) near_bnd(n_near_bnd) = c1
          exit
        end if
      end do

    end do
    if(run .eq. 1) allocate(near_bnd(n_near_bnd))
  end do
  print *, '# Total number of cells:         ', grid % n_cells
  print *, '# Number of near boundary cells: ', n_near_bnd

  !---------------------------------------------!
  !   Collect and sort data from inside cells   !
  !---------------------------------------------!
  n_bnd_proj = 0
  n_fr_cells = n_near_bnd / 20
  do cb = 1, n_near_bnd

    ! Get near boundary cell value
    c1 = near_bnd(cb)

    ! Get proper face
    if(grid % cells_n_nodes(c1) .eq. 4) fn = neu_tet
    if(grid % cells_n_nodes(c1) .eq. 5) fn = neu_pyr
    if(grid % cells_n_nodes(c1) .eq. 6) fn = neu_wed
    if(grid % cells_n_nodes(c1) .eq. 8) fn = neu_hex

    n_cell_faces = 6  ! assume it is a hexahedron
    if(grid % cells_n_nodes(c1) .eq. 4) n_cell_faces = 4
    if(grid % cells_n_nodes(c1) .eq. 5) n_cell_faces = 5
    if(grid % cells_n_nodes(c1) .eq. 6) n_cell_faces = 5

    ! Browse through all possible faces
    do j = 1, n_cell_faces
      n_match = 0

      n_face_nodes = 4         ! assume it is a quad
      if( fn(j, 4) < 0 ) then  ! nope, it was a triangle
        n_face_nodes = 3
      end if

      do n1 = 1, n_face_nodes
        if(fn(j, n1) > 0) then  ! if this node exists
          if( is_node_bnd( grid % cells_n(fn(j, n1), c1) ) ) then
            n_match = n_match + 1
          end if
        end if
      end do ! n1
      if(n_match .eq. n_face_nodes) then
        n_bnd_proj = n_bnd_proj + 1
        ai(1:4) = HUGE_INT
        ai(1:n_face_nodes) = grid % cells_n( fn(j,1:n_face_nodes), c1)
        call Sort_Mod_Int(ai)
        a1(n_bnd_proj) = ai(1)
        a2(n_bnd_proj) = ai(2)
        a3(n_bnd_proj) = ai(3)
        c (n_bnd_proj) = c1
        d (n_bnd_proj) = j
      end if ! n_match .eq. n_face_nodes
    end do ! j

  end do    ! c1(cb)
  call Sort_Mod_3_Int_Carry_2_Int(a1, a2, a3, c, d)

  !-----------------------------------------------!
  !   Collect and sort data from boundary cells   !
  !-----------------------------------------------!
  do c2 = -grid % n_bnd_cells, -1
    n_face_nodes = grid % cells_n_nodes(c2)
    ai(1:4) = HUGE_INT
    ai(1:n_face_nodes) = grid % cells_n(1:n_face_nodes, c2)
    call Sort_Mod_Int(ai)
    a1(-c2) = ai(1)
    a2(-c2) = ai(2)
    a3(-c2) = ai(3)
    b (-c2) = c2
  end do
  call Sort_Mod_3_Int_Carry_Int(a1, a2, a3, b)

  !---------------------------------!
  !   Finally, just match them up   !
  !---------------------------------!
  do cb = 1, n_bnd_proj
    grid % cells_bnd_color(d(cb), c(cb)) = grid % bnd_cond % color(b(cb))
  end do

  end subroutine
