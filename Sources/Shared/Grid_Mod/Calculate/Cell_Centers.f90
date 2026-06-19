!==============================================================================!
  subroutine Calculate_Cell_Centers(Grid)
!------------------------------------------------------------------------------!
!>  Calculate the cell centers from nodal coordinates.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type) :: Grid  !! grid under consideration
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, i_nod, n
!==============================================================================!

  !--------------------------------------------------------------------------!
  !   Compute them by adding node coordinates and dividing by their number   !
  !--------------------------------------------------------------------------!
  do c = 1, Grid % n_cells

    ! For concave cells, cell center was estimated in Convert_Mod/Create_Dual
    if(.not. Grid % concave(c)) then

      ! (Re)initialize cell coordinates
      Grid % xc(c) = 0.0
      Grid % yc(c) = 0.0
      Grid % zc(c) = 0.0

      do i_nod = 1, abs(Grid % cells_n_nodes(c))
        n = Grid % cells_n(i_nod, c)
        Grid % xc(c) = Grid % xc(c) + Grid % xn(n)  &
                     / real(abs(Grid % cells_n_nodes(c)))
        Grid % yc(c) = Grid % yc(c) + Grid % yn(n)  &
                     / real(abs(Grid % cells_n_nodes(c)))
        Grid % zc(c) = Grid % zc(c) + Grid % zn(n)  &
                     / real(abs(Grid % cells_n_nodes(c)))
      end do

    end if
  end do

  print *, '# Cell centers calculated !'

  end subroutine
