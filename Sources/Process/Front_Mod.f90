#include "../Shared/Assert.h90"
#include "../Shared/Browse.h90"

!==============================================================================!
  module Front_Mod
!------------------------------------------------------------------------------!
!   Module for Lagrangian particle tracking                                    !
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
  type Front_Type

    type(Grid_Type),  pointer :: pnt_grid  ! Grid for which it is defined
    type(Field_Type), pointer :: pnt_flow  ! flow field for which it is defined

    integer                      :: n_elems
    integer                      :: n_verts
    integer                      :: n_sides
    type(Vert_Type), allocatable :: Vert(:)
    type(Elem_Type), allocatable :: elem(:)
    type(Side_Type), allocatable :: side(:)

    ! Bounding nodes for each vertex (derives from Isoap usage)
    integer, allocatable :: b_node_1(:)
    integer, allocatable :: b_node_2(:)

    ! Is the mesh divided among processors (like in Front_Type)
    ! or is it shared among all of them (like in Surf_Type)
    logical :: mesh_divided

    ! Front to static Grid connectivity
    integer, allocatable :: elem_in_cell(:)
    integer, allocatable :: elems_at_face(:,:)

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
