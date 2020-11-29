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
  integer              :: fn(6,4), i_fac, i_nod, c1, c2, bc
  integer              :: n_match
  integer              :: n_bnd_proj   ! bnd. cells projected from cells
  integer              :: n_face_nodes ! number of nodes in a face
  integer              :: n_cell_faces ! number of faces in a cell
  logical, allocatable :: is_node_bnd(:)
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

  allocate(is_node_bnd(grid % n_nodes))

  n_bnd_proj = 0
  do bc = 1, grid % n_bnd_cond

    ! Re-initialize the logical array
    is_node_bnd(:) = .false.

    ! Mark boundary nodes in this section
    do c2 = -grid % n_bnd_cells, -1
      if( grid % bnd_cond % color(c2) .eq. bc ) then
        do i_nod = 1, grid % cells_n_nodes(c2)  ! 3 or 4
          is_node_bnd( grid % cells_n(i_nod, c2) ) = .true.
        end do
      end if
    end do

    do c1 = 1, grid % n_cells

      ! Fetch the faces from this cell
      if(grid % cells_n_nodes(c1) .eq. 4) fn = tet
      if(grid % cells_n_nodes(c1) .eq. 5) fn = pyr
      if(grid % cells_n_nodes(c1) .eq. 6) fn = wed
      if(grid % cells_n_nodes(c1) .eq. 8) fn = hex

      n_cell_faces = 6  ! assume it is a hexahedron
      if(grid % cells_n_nodes(c1) .eq. 4) n_cell_faces = 4
      if(grid % cells_n_nodes(c1) .eq. 5) n_cell_faces = 5
      if(grid % cells_n_nodes(c1) .eq. 6) n_cell_faces = 5

      ! Browse through all possible faces
      do i_fac = 1, n_cell_faces
        n_match = 0

        n_face_nodes = 4         ! assume it is a quad
        if( fn(i_fac, 4) < 0 ) then  ! nope, it was a triangle
          n_face_nodes = 3
        end if

        do i_nod = 1, n_face_nodes
          if(fn(i_fac, i_nod) > 0) then  ! if this node exists
            if( is_node_bnd( grid % cells_n(fn(i_fac, i_nod), c1) ) ) then
              n_match = n_match + 1
            end if
          end if
        end do ! i_nod

        ! This face is matching
        if(n_match .eq. n_face_nodes) then
          grid % cells_bnd_color(i_fac, c1) = bc
          n_bnd_proj = n_bnd_proj + 1
        end if ! n_match .eq. n_face_nodes

      end do ! i_fac
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
