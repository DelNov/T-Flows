!==============================================================================!
  subroutine Faces_Center(Grid, s, xf, yf, zf)
!------------------------------------------------------------------------------!
!>  Calculate a face's center from nodal coordinates.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type)    :: Grid        !! grid under consideration
  integer, intent(in) :: s           !! face number
  real,   intent(out) :: xf, yf, zf  !! face center cooridnates
!-----------------------------------[Locals]-----------------------------------!
  integer :: i_nod, n, nn
!==============================================================================!

  ! Nullify cooridinates
  xf = 0.0
  yf = 0.0
  zf = 0.0

  ! Copy face node coordinates to a local array for easier handling
  do i_nod = 1, Grid % faces_n_nodes(s)  ! local node counter
    n = Grid % faces_n(i_nod ,s)         ! global node number
    xf = xf + Grid % xn(n)
    yf = yf + Grid % yn(n)
    zf = zf + Grid % zn(n)
  end do

  ! Barycenter
  nn = Grid % faces_n_nodes(s)
  xf = xf / real(nn)
  yf = yf / real(nn)
  zf = zf / real(nn)

  end subroutine
