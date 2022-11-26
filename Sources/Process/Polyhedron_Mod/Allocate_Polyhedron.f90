!==============================================================================!
  subroutine Allocate_Polyhedron(Pol, n_faces, n_nodes)
!------------------------------------------------------------------------------!
!   Allocate memory for polyhedron                                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Polyhedron_Type) :: Pol
  integer, intent(in)    :: n_faces
  integer, intent(in)    :: n_nodes
!==============================================================================!

  allocate(Pol % faces_n_nodes(n_faces))
  allocate(Pol % faces_n      (n_faces, n_nodes))
  allocate(Pol % nodes_xyz    (n_nodes, 3))
  allocate(Pol % phi          (n_nodes))
  allocate(Pol % phi_int      (n_nodes))
  allocate(Pol % global_node  (n_nodes))

  end subroutine
