!==============================================================================!
  subroutine Grid_Mod_Calculate_Cell_Centers(grid)
!------------------------------------------------------------------------------!
!   Calculate the cell centers from nodal coordinates.                         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, i_nod, n
!==============================================================================!

  ! (Re)initialize cell coordinates
  grid % xc(:) = 0.0
  grid % yc(:) = 0.0
  grid % zc(:) = 0.0

  ! Compute them by adding node coordinates and dividing by their number
  do c = 1, grid % n_cells
    do i_nod = 1, abs(grid % cells_n_nodes(c))
      n = grid % cells_n(i_nod, c)
      grid % xc(c) = grid % xc(c) + grid % xn(n)  &
                   / real(abs(grid % cells_n_nodes(c)))
      grid % yc(c) = grid % yc(c) + grid % yn(n)  &
                   / real(abs(grid % cells_n_nodes(c)))
      grid % zc(c) = grid % zc(c) + grid % zn(n)  &
                   / real(abs(grid % cells_n_nodes(c)))
    end do
  end do

  print *, '# Cell centers calculated !'

  end subroutine
