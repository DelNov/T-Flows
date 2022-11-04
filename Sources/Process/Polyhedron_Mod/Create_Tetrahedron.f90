!==============================================================================!
  subroutine Create_Tetrahedron(Pol)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Polyhedron_Type) :: Pol
!==============================================================================!

  Pol % n_faces = 4
  Pol % n_nodes = 4
  Pol % faces_n_nodes(1) = 3
  Pol % faces_n(1,1) = 1
  Pol % faces_n(1,2) = 3
  Pol % faces_n(1,3) = 2
  Pol % faces_n_nodes(2) = 3
  Pol % faces_n(2,1) = 3
  Pol % faces_n(2,2) = 1
  Pol % faces_n(2,3) = 4
  Pol % faces_n_nodes(3) = 3
  Pol % faces_n(3,1) = 2
  Pol % faces_n(3,2) = 3
  Pol % faces_n(3,3) = 4
  Pol % faces_n_nodes(4) = 3
  Pol % faces_n(4,1) = 1
  Pol % faces_n(4,2) = 2
  Pol % faces_n(4,3) = 4
  Pol % nodes_xyz(1,1) = 0.0
  Pol % nodes_xyz(1,2) = 0.0
  Pol % nodes_xyz(1,3) = 0.0
  Pol % nodes_xyz(2,1) = 0.91
  Pol % nodes_xyz(2,2) = 0.24
  Pol % nodes_xyz(2,3) = 1.0
  Pol % nodes_xyz(3,1) = 0.72
  Pol % nodes_xyz(3,2) = 0.16
  Pol % nodes_xyz(3,3) = 0.07
  Pol % nodes_xyz(4,1) = 1.0
  Pol % nodes_xyz(4,2) = 1.0
  Pol % nodes_xyz(4,3) = 1.0

  end
