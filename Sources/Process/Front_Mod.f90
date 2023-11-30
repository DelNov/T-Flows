#include "../Shared/Assert.h90"
#include "../Shared/Browse.h90"

!==============================================================================!
  module Front_Mod
!------------------------------------------------------------------------------!
!>  This module is designed for generating and managing a front, typically an
!>  interface between two phases in a VOF context. This module constructs the
!>  front using polygons spanned across the numerical mesh, thus it can be said
!>  that it is "mesh-dependent". It places polygon vertices at cell edges,
!>  aligning them these with specific field variable values, typically the 0.5
!>  iso-surface in a VOF context, where values range from 0 to 1. Through
!>  Front_Type which is defined in this module, it introduces a range of
!>  functionalities for managing and saving a front in .vtu format for
!>  post-processing.  The principal functionality of this module, to extract
!>  surface polygons based on field functions, is based on Isoap library
!>  (https://data.mendeley.com/datasets/4rcf98s74c) and through Isoap_Mod.
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Field_Mod
  use Vert_Mod
  use Elem_Mod
  use Side_Mod
  use Isoap_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !----------------!
  !   Front type   !
  !----------------!
  !> This type incorporates various elements, vertices, and sides to form a
  !> mesh that aligns with specific field values. It also holds functionalities
  !> for creating and manipulating front and saving it in .vtu format.
  type Front_Type

    type(Grid_Type),  pointer :: pnt_grid  !! grid for which it is defined
    type(Field_Type), pointer :: pnt_flow  !! flow field for which it is defined

    integer                      :: n_elems  !! number of elements
    integer                      :: n_verts  !! number of vertices
    integer                      :: n_sides  !! number of sides
    type(Vert_Type), allocatable :: Vert(:)  !! list of vertices
    type(Elem_Type), allocatable :: Elem(:)  !! list of elements
    type(Side_Type), allocatable :: side(:)  !! list of sides

    ! Bounding nodes for each vertex (derives from Isoap usage)
    integer, allocatable :: b_node_1(:)
    integer, allocatable :: b_node_2(:)

    ! Is the mesh divided among processors (like in Front_Type)
    ! or is it shared among all of them (like in Surf_Type)
    logical :: mesh_divided

    ! Front to static Grid connectivity
    integer, allocatable :: elem_in_cell(:)

    ! Face-based interserction with surface
    logical, allocatable :: intersects_face(:)
    real,    allocatable :: xs(:), ys(:), zs(:)

    contains

      procedure :: Allocate_Front
      procedure :: Calculate_Element_Centroids
      procedure :: Calculate_Element_Normals
      procedure :: Check_Elements
      procedure :: Clean_Front
      procedure :: Compress_Front_Vertices
      procedure :: Compute_Distance_Function_And_Vof
      procedure :: Find_Front_Elements
      procedure :: Find_Sides
      procedure :: Find_Vertex_Elements
      procedure :: Initialize_Front
      procedure :: Mark_Cells_And_Faces
      procedure :: Place_Front_At_Value
      procedure :: Print_Front_Statistics
      procedure :: Save_Debug_Front_Vtu

  end type

  contains

#   include "Front_Mod/Allocate_Front.f90"
#   include "Front_Mod/Calculate_Element_Centroids.f90"
#   include "Front_Mod/Calculate_Element_Normals.f90"
#   include "Front_Mod/Check_Elements.f90"
#   include "Front_Mod/Clean.f90"
#   include "Front_Mod/Compress_Front_Vertices.f90"
#   include "Front_Mod/Compute_Distance_Function_And_Vof.f90"
#   include "Front_Mod/Find_Front_Elements.f90"
#   include "Front_Mod/Find_Sides.f90"
#   include "Front_Mod/Find_Vertex_Elements.f90"
#   include "Front_Mod/Initialize_Front.f90"
#   include "Front_Mod/Mark_Cells_And_Faces.f90"
#   include "Front_Mod/Place_Front_At_Value.f90"
#   include "Front_Mod/Print_Front_Statistics.f90"
#   include "Front_Mod/Save_Debug_Front_Vtu.f90"

  end module
