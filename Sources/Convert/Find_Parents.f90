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
  integer :: fn(6,4), j, n1, n2, c1, c2
  integer :: n_match, number_of_found_parents, n_face_nodes, n_cells_fraction
!==============================================================================!

  print *, '#================================================='
  print *, '# Parent information not given in the input file!'
  print *, '# Looking for parents. This may take a few minutes'
  print *, '#-------------------------------------------------'

  number_of_found_parents = 0
  n_cells_fraction = grid % n_cells / 20  ! each 5%
  do c1 = 1, grid % n_cells

    ! Fetch the right shape
    if(grid % cells_n_nodes(c1) .eq. 4) fn = neu_tet
    if(grid % cells_n_nodes(c1) .eq. 5) fn = neu_pyr
    if(grid % cells_n_nodes(c1) .eq. 6) fn = neu_wed
    if(grid % cells_n_nodes(c1) .eq. 8) fn = neu_hex

    ! Print some statistics on the screen
    if(mod( c1, n_cells_fraction ) .eq. 0) then
      print '(a2, f5.0, a14)',              &
        ' #',                               &
        (100. * c1 / (1.0*(grid % n_cells))),  &
        ' % complete...'
    end if ! each 5%

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
          number_of_found_parents = number_of_found_parents + 1
        end if ! n_match .eq. n_face_nodes
      end do ! j

    end do
  end do

  print *, '# Number of found parents = ', number_of_found_parents

  end subroutine
