!==============================================================================!
  subroutine Cells_Center(Grid, c, xc, yc, zc)
!------------------------------------------------------------------------------!
!>  Calculate a cell's center from nodal coordinates.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type)    :: Grid        !! grid under consideration
  integer, intent(in) :: c           !! face number
  real,   intent(out) :: xc, yc, zc  !! face center cooridnates
!-----------------------------------[Locals]-----------------------------------!
  integer :: i_nod, n, nn
!==============================================================================!

  ! Nullify cooridinates
  xc = 0.0
  yc = 0.0
  zc = 0.0

  ! Browse throug cell's nodes
  do i_nod = 1, Grid % cells_n_nodes(c)  ! local node counter
    n = Grid % cells_n(i_nod, c)         ! global node number
    xc = xc + Grid % xn(n)
    yc = yc + Grid % yn(n)
    zc = zc + Grid % zn(n)
  end do

  ! Barycenter
  nn = Grid % cells_n_nodes(c)
  xc = xc / real(nn)
  yc = yc / real(nn)
  zc = zc / real(nn)

  end subroutine
