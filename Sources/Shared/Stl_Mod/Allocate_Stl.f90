!==============================================================================!
  subroutine Allocate_Stl(Stl, n_facets)
!------------------------------------------------------------------------------!
!>  This subroutine initializes an Stl object, part of the Stl_Mod, for the
!>  storage of geometry in the STL format. It allocates memory for storing
!>  facet details such as normals and coordinates, setting the stage for
!>  further processing of STL files.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Stl_Type)     :: Stl       !! parent Stl_Type object
  integer, intent(in) :: n_facets  !! number of facets in STL
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
