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
  integer :: fn(6,4), j, n1, n2, c1, c2, n_match, number_of_found_parents
!==============================================================================!

  print *, '# Parent information not given in CGNS file!'
  print *, '# Looking for parents. This may take a few minutes'

  number_of_found_parents = 0
  do c1 = 1, grid % n_cells
    if(grid % cells_n_nodes(c1) .eq. 4) fn = neu_tet
    if(grid % cells_n_nodes(c1) .eq. 5) fn = neu_pyr
    if(grid % cells_n_nodes(c1) .eq. 6) fn = neu_wed
    if(grid % cells_n_nodes(c1) .eq. 8) fn = neu_hex
    do c2 = -grid % n_bnd_cells, -1

      !------------------------------!
      !   Number of matching nodes   !
      !------------------------------!
      do j = 1, 6     ! six is maximum number of faces (hexahedron)
        n_match = 0
        do n1 = 1, 4  ! four is maximum number of nodes (quad)
          do n2 = 1, grid % cells_n_nodes(c2)
            if(fn(j, n1) > 0) then  ! if this node exists
              if(grid % cells_n(fn(j, n1), c1) .eq.  &
                 grid % cells_n(n2,        c2)) then
                n_match = n_match + 1
              end if
            end if
          end do
        end do
        if(n_match >= 3) then
          grid % cells_bnd_color(j, c1) = grid % bnd_cond % color(c2)
          number_of_found_parents = number_of_found_parents + 1
        end if
      end do

    end do
  end do

  end subroutine
