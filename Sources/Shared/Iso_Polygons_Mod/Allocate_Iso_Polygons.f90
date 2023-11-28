!==============================================================================!
  subroutine Allocate_Iso_Polygons(Iso, n_faces, n_nodes)
!------------------------------------------------------------------------------!
!>  This subroutine allocates memory for various arrays within the Iso object
!>  (an object of the type Iso_Polygons_Type).  After allocation, it sets the
!>  allocated flag of Iso to .true.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Iso_Polygons_Type) :: Iso      !! parent class
  integer, intent(in)      :: n_faces  !! number of faces
  integer, intent(in)      :: n_nodes  !! number of nodes
!==============================================================================!

  allocate(Iso % polys_v      (n_faces, n_nodes));
  allocate(Iso % face_index   (n_nodes))
  allocate(Iso % polys_n_verts(n_faces))
  allocate(Iso % verts_xyz    (n_nodes, 3))
  allocate(Iso % b_node_1     (n_nodes))
  allocate(Iso % b_node_2     (n_nodes))
  Iso % allocated = .true.

  end subroutine
