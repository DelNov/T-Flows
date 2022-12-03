!==============================================================================!
  subroutine Create_Pentapyramid(Pol)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Polyhedron_Type) :: Pol
!==============================================================================!

  Pol % n_faces = 6
  Pol % n_nodes = 6

  Pol % nodes_xyz(1,1) = 0.04
  Pol % nodes_xyz(1,2) = 0.77
  Pol % nodes_xyz(1,3) = 0.0
  Pol % nodes_xyz(2,1) = 0.0
  Pol % nodes_xyz(2,2) = 0.0
  Pol % nodes_xyz(2,3) = 0.0
  Pol % nodes_xyz(3,1) = 0.49
  Pol % nodes_xyz(3,2) = 0.22
  Pol % nodes_xyz(3,3) = 0.0
  Pol % nodes_xyz(4,1) = 1.0
  Pol % nodes_xyz(4,2) = 0.13
  Pol % nodes_xyz(4,3) = 0.0
  Pol % nodes_xyz(5,1) = 0.16
  Pol % nodes_xyz(5,2) = 1.0
  Pol % nodes_xyz(5,3) = 0.0
  Pol % nodes_xyz(6,1) = 0.1
  Pol % nodes_xyz(6,2) = 0.5
  Pol % nodes_xyz(6,3) = 1.0

  ! is = 1
  Pol % faces_n_nodes(1) = 5
  Pol % faces_n(1,1) = 1
  Pol % faces_n(1,2) = 5
  Pol % faces_n(1,3) = 4
  Pol % faces_n(1,4) = 3
  Pol % faces_n(1,5) = 2
  ! is = 2
  Pol % faces_n_nodes(2) = 3
  Pol % faces_n(2,1) = 1
  Pol % faces_n(2,2) = 2
  Pol % faces_n(2,3) = 6
  ! is = 3
  Pol % faces_n_nodes(3) = 3
  Pol % faces_n(3,1) = 2
  Pol % faces_n(3,2) = 3
  Pol % faces_n(3,3) = 6
  ! is = 4
  Pol % faces_n_nodes(4) = 3
  Pol % faces_n(4,1) = 3
  Pol % faces_n(4,2) = 4
  Pol % faces_n(4,3) = 6
  ! is = 5
  Pol % faces_n_nodes(5) = 3
  Pol % faces_n(5,1) = 4
  Pol % faces_n(5,2) = 5
  Pol % faces_n(5,3) = 6
  ! is = 6
  Pol % faces_n_nodes(6) = 3
  Pol % faces_n(6,1) = 5
  Pol % faces_n(6,2) = 1
  Pol % faces_n(6,3) = 6

  end subroutine
