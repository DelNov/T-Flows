!==============================================================================!
  subroutine Allocate_Iso_Polygons(Iso, n_faces, n_nodes)
!------------------------------------------------------------------------------!
!   Allocate memory for iso-polygons                                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Iso_Polygons_Type) :: Iso
  integer, intent(in)      :: n_faces
  integer, intent(in)      :: n_nodes
!==============================================================================!

  allocate(Iso % polys_v      (n_faces, n_nodes));
  allocate(Iso % face_index   (n_nodes))
  allocate(Iso % polys_n_verts(n_faces))
  allocate(Iso % verts_xyz    (n_nodes, 3))
  allocate(Iso % b_node_1     (n_nodes))
  allocate(Iso % b_node_2     (n_nodes))

  end subroutine
