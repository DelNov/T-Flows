!==============================================================================!
  pure subroutine Calculate_Cell_Centroid(Pol, xc, yc, zc)
!------------------------------------------------------------------------------!
!>  This subroutine is designed to calculate polyhedral cell's centroid.
!>  This information is not stored in the Polyhedron_Type object.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Polyhedron_Type), intent(in)  :: Pol         !! parent class
  real,                   intent(out) :: xc, yc, zc  !! centroid's coordinate
!-----------------------------------[Locals]-----------------------------------!
  integer :: i
!==============================================================================!

  xc = 0.0
  yc = 0.0
  zc = 0.0
  do i = 1, Pol % n_nodes
    xc = xc + Pol % nodes_xyz(i, 1)
    yc = yc + Pol % nodes_xyz(i, 2)
    zc = zc + Pol % nodes_xyz(i, 3)
  end do
  xc = xc / real(Pol % n_nodes)
  yc = yc / real(Pol % n_nodes)
  zc = zc / real(Pol % n_nodes)

  end subroutine
