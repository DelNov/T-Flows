!==============================================================================!
  subroutine Grid_Mod_Calculate_Cell_Centers(grid)
!------------------------------------------------------------------------------!
!   Calculate the cell centers from nodal coordinates.                         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, n
!==============================================================================!

  ! (Re)initialize cell coordinates
  grid % xc(:) = 0.0
  grid % yc(:) = 0.0
  grid % zc(:) = 0.0

  ! Compute them by adding node coordinates and dividing by their number
  do c = 1, grid % n_cells
    do n = 1, grid % cells_n_nodes(c)
      grid % xc(c) = grid % xc(c) + grid % xn(grid % cells_n(n,c))  &
                   / (1.0*grid % cells_n_nodes(c))
      grid % yc(c) = grid % yc(c) + grid % yn(grid % cells_n(n,c))  &
                   / (1.0*grid % cells_n_nodes(c))
      grid % zc(c) = grid % zc(c) + grid % zn(grid % cells_n(n,c))  &
                   / (1.0*grid % cells_n_nodes(c))
    end do
  end do

  print *, '# Cell centers calculated !'

  end subroutine
