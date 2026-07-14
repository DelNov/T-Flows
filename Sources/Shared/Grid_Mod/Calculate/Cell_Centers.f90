!==============================================================================!
  subroutine Calculate_Cell_Centers(Grid)
!------------------------------------------------------------------------------!
!>  Calculate the cell centers from nodal coordinates.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type) :: Grid  !! grid under consideration
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, i_nod, n, nn
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

      nn = abs(Grid % cells_n_nodes(c))
      Assert(nn .gt. 0)

      do i_nod = 1, nn
        n = Grid % cells_n(i_nod, c)
        Grid % xc(c) = Grid % xc(c) + Grid % xn(n)
        Grid % yc(c) = Grid % yc(c) + Grid % yn(n)
        Grid % zc(c) = Grid % zc(c) + Grid % zn(n)
      end do

      ! Barycenter
      Grid % xc(c) = Grid % xc(c) / real(nn)
      Grid % yc(c) = Grid % yc(c) / real(nn)
      Grid % zc(c) = Grid % zc(c) / real(nn)

    end if

  end do

  print *, '# Cell centers calculated !'

  end subroutine
