!==============================================================================!
  subroutine Find_Parents(grid)
!------------------------------------------------------------------------------!
!   Looks boundary cells' parents for meshes in which they are not given       !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod, only: Grid_Type
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!------------------------------------------------------------------------------!
  include 'Cell_Numbering_Neu.f90'
!-----------------------------------[Locals]-----------------------------------!
  integer              :: fn(6,4), j, n1, n2, c1, c2, cb, run
  integer              :: n_match, n_found_parents, n_cells_fraction
  integer              :: n_bnd_nodes  ! near boundary nodes
  integer              :: n_bnd_proj   ! bnd. cells projected from cells
  integer              :: n_near_bnd   ! number of near boundary cells
  integer              :: n_face_nodes ! number of nodes in a face
  integer              :: n_fr_cells   ! for statistics
  integer, allocatable :: near_bnd(:)  ! near boundary cells
  logical, allocatable :: is_node_bnd(:)
  logical, allocatable :: is_cell_marked(:)
!==============================================================================!

  print *, '#================================================='
  print *, '# Parent information not given in CGNS file!'
  print *, '# Looking for parents. This may take a few minutes'
  print *, '#-------------------------------------------------'

  allocate(is_node_bnd   (grid % n_nodes))
  allocate(is_cell_marked(grid % n_cells))
  is_node_bnd   (:) = .false.
  is_cell_marked(:) = .false.

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

  n_bnd_proj = 0
  n_fr_cells = n_near_bnd / 20
  do cb = 1, n_near_bnd

    ! Print some statistics on the screen
    if(mod( cb, n_fr_cells ) .eq. 0) then
      print '(a2, f5.0, a14)',              &
        ' #',                               &
        (100. * cb / (1.0*(n_near_bnd))),  &
        ' % complete...'
    end if ! each 5%

    ! Get near boundary cell value
    c1 = near_bnd(cb)

    ! Get proper face
    if(grid % cells_n_nodes(c1) .eq. 4) fn = neu_tet
    if(grid % cells_n_nodes(c1) .eq. 5) fn = neu_pyr
    if(grid % cells_n_nodes(c1) .eq. 6) fn = neu_wed
    if(grid % cells_n_nodes(c1) .eq. 8) fn = neu_hex

    do c2 = -grid % n_bnd_cells, -1

      !------------------------------!
      !   Number of matching nodes   !
      !------------------------------!
      do j = 1, 6  ! from 1 to 6th face (6 = hexahedron)
        n_match = 0

        n_face_nodes = 4         ! assume it is a quad
        if( fn(j, 4) < 0 ) then  ! nope, it was a triangle
          n_face_nodes = 3
        end if

        do n1 = 1, n_face_nodes
          do n2 = 1, grid % cells_n_nodes(c2)  ! from 1 to 3/4th node in face
            if(fn(j, n1) > 0) then  ! if this node exists
              if(grid % cells_n(fn(j, n1), c1) .eq.  &
                 grid % cells_n(n2,        c2)) then
                n_match = n_match + 1
              end if
            end if
          end do ! n2
        end do ! n1
        if(n_match .eq. n_face_nodes) then
          grid % cells_bnd_color(j, c1) = grid % bnd_cond % color(c2)
          n_bnd_proj = n_bnd_proj + 1
        end if ! n_match .eq. n_face_nodes
      end do ! j

    end do  ! c2
  end do    ! c1(cb)

  print *, '# Number of boundary cells:           ', grid % n_bnd_cells
  print *, '# Number of projected boundary cells: ', n_bnd_proj

  end subroutine
