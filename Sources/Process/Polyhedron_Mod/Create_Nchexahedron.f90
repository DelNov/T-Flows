!==============================================================================!
  subroutine Create_Nchexahedron(Pol)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Polyhedron_Type) :: Pol
!-----------------------------------[Locals]-----------------------------------!
  real :: d0, d1
!==============================================================================!

  d0 = 0.0
  d1 = 1.0

  Pol % n_faces = 6
  Pol % n_nodes = 8
  Pol % faces_n_nodes(1) = 4
  Pol % faces_n(1,1) = 1
  Pol % faces_n(1,2) = 2
  Pol % faces_n(1,3) = 3
  Pol % faces_n(1,4) = 4
  Pol % faces_n_nodes(2) = 4
  Pol % faces_n(2,1) = 2
  Pol % faces_n(2,2) = 1
  Pol % faces_n(2,3) = 5
  Pol % faces_n(2,4) = 6
  Pol % faces_n_nodes(3) = 4
  Pol % faces_n(3,1) = 3
  Pol % faces_n(3,2) = 2
  Pol % faces_n(3,3) = 6
  Pol % faces_n(3,4) = 7
  Pol % faces_n_nodes(4) = 4
  Pol % faces_n(4,1) = 4
  Pol % faces_n(4,2) = 3
  Pol % faces_n(4,3) = 7
  Pol % faces_n(4,4) = 8
  Pol % faces_n_nodes(5) = 4
  Pol % faces_n(5,1) = 1
  Pol % faces_n(5,2) = 4
  Pol % faces_n(5,3) = 8
  Pol % faces_n(5,4) = 5
  Pol % faces_n_nodes(6) = 4
  Pol % faces_n(6,1) = 6
  Pol % faces_n(6,2) = 5
  Pol % faces_n(6,3) = 8
  Pol % faces_n(6,4) = 7
  Pol % nodes_xyz(1,1) = 0.5
  Pol % nodes_xyz(1,2) = 0.75
  Pol % nodes_xyz(1,3) = d1
  Pol % nodes_xyz(2,1) = 0.5
  Pol % nodes_xyz(2,2) = 0.75
  Pol % nodes_xyz(2,3) = d0
  Pol % nodes_xyz(3,1) = d1
  Pol % nodes_xyz(3,2) = d1
  Pol % nodes_xyz(3,3) = d0
  Pol % nodes_xyz(4,1) = d1
  Pol % nodes_xyz(4,2) = d1
  Pol % nodes_xyz(4,3) = d1
  Pol % nodes_xyz(5,1) = d0
  Pol % nodes_xyz(5,2) = d0
  Pol % nodes_xyz(5,3) = d1
  Pol % nodes_xyz(6,1) = d0
  Pol % nodes_xyz(6,2) = d0
  Pol % nodes_xyz(6,3) = d0
  Pol % nodes_xyz(7,1) = d0
  Pol % nodes_xyz(7,2) = d1
  Pol % nodes_xyz(7,3) = d0
  Pol % nodes_xyz(8,1) = d0
  Pol % nodes_xyz(8,2) = d1
  Pol % nodes_xyz(8,3) = d1

  end
