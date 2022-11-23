!==============================================================================!
  subroutine Calculate_Cell_Centroid(Polyhedron, cell, xc, yc, zc)
!------------------------------------------------------------------------------!
!   Calculate cell's centroid, this information is not stored in Polyhedron    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Polyhedron_Type) :: Polyhedron
  integer, intent(in)    :: cell
  real,    intent(out)   :: xc, yc, zc
!-----------------------------------[Locals]-----------------------------------!
  integer :: i
!==============================================================================!

  xc = 0.0
  yc = 0.0
  zc = 0.0
  do i = 1, Polyhedron % n_nodes
    xc = xc + Polyhedron % nodes_xyz(i, 1)
    yc = yc + Polyhedron % nodes_xyz(i, 2)
    zc = zc + Polyhedron % nodes_xyz(i, 3)
  end do
  xc = xc / real(Polyhedron % n_nodes)
  yc = yc / real(Polyhedron % n_nodes)
  zc = zc / real(Polyhedron % n_nodes)

  end subroutine
