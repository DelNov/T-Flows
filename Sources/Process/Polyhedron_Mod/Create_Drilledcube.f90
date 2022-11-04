!==============================================================================!
  subroutine Create_Drilledcube(Pol)
!------------------------------------------------------------------------------!
!   Non-simply connected polyhedron                                            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Polyhedron_Type) :: Pol
!-----------------------------------[Locals]-----------------------------------!
  real :: d0, d02, d1, d12
!==============================================================================!

  call Pol % Create_Cube()

  d0 = 0.0
  d1 = 1.0

  Pol % faces_n_nodes(Pol % n_faces+1) = 4
  Pol % faces_n(Pol % n_faces+1,4) = Pol % n_nodes+1
  Pol % faces_n(Pol % n_faces+1,3) = Pol % n_nodes+2
  Pol % faces_n(Pol % n_faces+1,2) = Pol % n_nodes+3
  Pol % faces_n(Pol % n_faces+1,1) = Pol % n_nodes+4
  Pol % faces_n_nodes(Pol % n_faces+2) = 4
  Pol % faces_n(Pol % n_faces+2,4) = Pol % n_nodes+2
  Pol % faces_n(Pol % n_faces+2,3) = Pol % n_nodes+1
  Pol % faces_n(Pol % n_faces+2,2) = Pol % n_nodes+5
  Pol % faces_n(Pol % n_faces+2,1) = Pol % n_nodes+6
  Pol % faces_n_nodes(Pol % n_faces+3) = 4
  Pol % faces_n(Pol % n_faces+3,4) = Pol % n_nodes+3
  Pol % faces_n(Pol % n_faces+3,3) = Pol % n_nodes+2
  Pol % faces_n(Pol % n_faces+3,2) = Pol % n_nodes+6
  Pol % faces_n(Pol % n_faces+3,1) = Pol % n_nodes+7
  Pol % faces_n_nodes(Pol % n_faces+4) = 4
  Pol % faces_n(Pol % n_faces+4,4) = Pol % n_nodes+4
  Pol % faces_n(Pol % n_faces+4,3) = Pol % n_nodes+3
  Pol % faces_n(Pol % n_faces+4,2) = Pol % n_nodes+7
  Pol % faces_n(Pol % n_faces+4,1) = Pol % n_nodes+8
  Pol % faces_n_nodes(Pol % n_faces+5) = 4
  Pol % faces_n(Pol % n_faces+5,4) = Pol % n_nodes+1
  Pol % faces_n(Pol % n_faces+5,3) = Pol % n_nodes+4
  Pol % faces_n(Pol % n_faces+5,2) = Pol % n_nodes+8
  Pol % faces_n(Pol % n_faces+5,1) = Pol % n_nodes+5
  Pol % faces_n_nodes(Pol % n_faces+6) = 4
  Pol % faces_n(Pol % n_faces+6,4) = Pol % n_nodes+6
  Pol % faces_n(Pol % n_faces+6,3) = Pol % n_nodes+5
  Pol % faces_n(Pol % n_faces+6,2) = Pol % n_nodes+8
  Pol % faces_n(Pol % n_faces+6,1) = Pol % n_nodes+7

  d02 = 0.25
  d12 = 0.75

  Pol % nodes_xyz(Pol % n_nodes+1,1) = d12
  Pol % nodes_xyz(Pol % n_nodes+1,2) = d0
  Pol % nodes_xyz(Pol % n_nodes+1,3) = d12
  Pol % nodes_xyz(Pol % n_nodes+2,1) = d12
  Pol % nodes_xyz(Pol % n_nodes+2,2) = d0
  Pol % nodes_xyz(Pol % n_nodes+2,3) = d02
  Pol % nodes_xyz(Pol % n_nodes+3,1) = d12
  Pol % nodes_xyz(Pol % n_nodes+3,2) = d1
  Pol % nodes_xyz(Pol % n_nodes+3,3) = d02
  Pol % nodes_xyz(Pol % n_nodes+4,1) = d12
  Pol % nodes_xyz(Pol % n_nodes+4,2) = d1
  Pol % nodes_xyz(Pol % n_nodes+4,3) = d12
  Pol % nodes_xyz(Pol % n_nodes+5,1) = d02
  Pol % nodes_xyz(Pol % n_nodes+5,2) = d0
  Pol % nodes_xyz(Pol % n_nodes+5,3) = d12
  Pol % nodes_xyz(Pol % n_nodes+6,1) = d02
  Pol % nodes_xyz(Pol % n_nodes+6,2) = d0
  Pol % nodes_xyz(Pol % n_nodes+6,3) = d02
  Pol % nodes_xyz(Pol % n_nodes+7,1) = d02
  Pol % nodes_xyz(Pol % n_nodes+7,2) = d1
  Pol % nodes_xyz(Pol % n_nodes+7,3) = d02
  Pol % nodes_xyz(Pol % n_nodes+8,1) = d02
  Pol % nodes_xyz(Pol % n_nodes+8,2) = d1
  Pol % nodes_xyz(Pol % n_nodes+8,3) = d12

  Pol % n_faces = Pol % n_faces+6
  Pol % n_nodes = Pol % n_nodes+8

  end
