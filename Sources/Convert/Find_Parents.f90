!==============================================================================!
  subroutine Find_Parents(grid)
!------------------------------------------------------------------------------!
!   Looks boundary cells' parents for meshes in which they are not given       !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod,  only: Grid_Type
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!------------------------------------------------------------------------------!
  include 'Cell_Numbering_Neu.f90'
!-----------------------------------[Locals]-----------------------------------!
  integer              :: fn(6,4), j, n1, n2, c1, c2, bc
  integer              :: n_match
  integer              :: n_bnd_proj   ! bnd. cells projected from cells
  integer              :: n_face_nodes ! number of nodes in a face
  integer              :: n_cell_faces ! number of faces in a cell
  logical, allocatable :: is_node_bnd(:)
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

  allocate(is_node_bnd(grid % n_nodes))

  n_bnd_proj = 0
  do bc = 1, grid % n_bnd_cond

    ! Re-initialize the logical array
    is_node_bnd(:) = .false.

    ! Mark boundary nodes in this section
    do c2 = -grid % n_bnd_cells, -1
      if( grid % bnd_cond % color(c2) .eq. bc ) then
        do n2 = 1, grid % cells_n_nodes(c2)  ! 3 or 4
          is_node_bnd( grid % cells_n(n2, c2) ) = .true.
        end do
      end if
    end do

    do c1 = 1, grid % n_cells

      ! Fetch the faces from this cell
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

        ! This face is matching
        if(n_match .eq. n_face_nodes) then
          grid % cells_bnd_color(j, c1) = bc
          n_bnd_proj = n_bnd_proj + 1
        end if ! n_match .eq. n_face_nodes

      end do ! j
    end do
  end do

  !-------------------------------------------!
  !   Update number of boundary cells         !
  !   (This is probably not needed because    !
  !   it is re-calculated in Grid_Topology)   !
  !-------------------------------------------!
  if(n_bnd_proj > grid % n_bnd_cells) then
    grid % n_bnd_cells = n_bnd_proj
  end if

  end subroutine
