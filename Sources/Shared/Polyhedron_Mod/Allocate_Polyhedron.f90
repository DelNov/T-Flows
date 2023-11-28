!==============================================================================!
  subroutine Allocate_Polyhedron(Pol, n_faces, n_nodes)
!------------------------------------------------------------------------------!
!>  This subroutine allocates memory for various arrays within the Pol object
!>  (object of Polyhedron_Type). After allocation, it sets the allocated flag
!>  of Pol to .true.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Polyhedron_Type), intent(out) :: Pol      !! parent class
  integer,                intent(in)  :: n_faces  !! number of faces
  integer,                intent(in)  :: n_nodes  !! number of nodes
!==============================================================================!

  allocate(Pol % faces_n_nodes(n_faces))
  allocate(Pol % faces_n      (n_faces, n_nodes))
  allocate(Pol % nodes_xyz    (n_nodes, 3))
  allocate(Pol % phi          (n_nodes))
  allocate(Pol % global_node  (n_nodes))
  Pol % allocated = .true.

  end subroutine
