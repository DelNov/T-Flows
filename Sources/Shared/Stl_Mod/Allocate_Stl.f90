!==============================================================================!
  subroutine Allocate_Stl(Stl, n_facets)
!------------------------------------------------------------------------------!
!   Creates an Stl object from a file                                          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Stl_Type)     :: Stl
  integer, intent(in) :: n_facets
!==============================================================================!

  Stl % n_bnd_cells =  n_facets
  Stl % n_nodes     =  Stl % n_bnd_cells * 3
  Stl % n_cells     = -1
  Stl % n_faces     =  0
  Stl % n_edges     =  0
  allocate(Stl % xn(Stl % n_nodes));          Stl % xn(:) = 0.0
  allocate(Stl % yn(Stl % n_nodes));          Stl % yn(:) = 0.0
  allocate(Stl % zn(Stl % n_nodes));          Stl % zn(:) = 0.0
  allocate(Stl % nx(-Stl % n_bnd_cells:-1));  Stl % nx(:) = 0.0
  allocate(Stl % ny(-Stl % n_bnd_cells:-1));  Stl % ny(:) = 0.0
  allocate(Stl % nz(-Stl % n_bnd_cells:-1));  Stl % nz(:) = 0.0
  allocate(Stl % cells_n_nodes(-Stl % n_bnd_cells:-1))
  allocate(Stl % cells_n   (3, -Stl % n_bnd_cells:-1))
  Stl % cells_n_nodes(:) = 3
  Stl % cells_n(:,:)     = 0

  end subroutine
