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

  ! For concave cells, cell center was estimated in Convert_Mod/Create_Dual
  if(.not. Grid % concave(c)) then

    nn = abs(Grid % cells_n_nodes(c))
    Assert(nn .gt. 0)

    ! Browse throug cell's nodes
    do i_nod = 1, nn                ! local node counter
      n = Grid % cells_n(i_nod, c)  ! global node number
      xc = xc + Grid % xn(n)
      yc = yc + Grid % yn(n)
      zc = zc + Grid % zn(n)
    end do

    ! Barycenter
    xc = xc / real(nn)
    yc = yc / real(nn)
    zc = zc / real(nn)

  ! Cell is concave, take the already computed cell center
  else

    xc = Grid % xc(c)
    yc = Grid % yc(c)
    zc = Grid % zc(c)

  end if


  end subroutine
