!==============================================================================!
  module Front_Mod
!------------------------------------------------------------------------------!
!   Module for Lagrangian particle tracking                                    !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Profiler_Mod
  use Field_Mod
  use Vert_Mod
  use Elem_Mod
  use Side_Mod
  use Isoap_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !--------------------------------!
  !   A few important parameters   !
  !--------------------------------!
  integer, parameter   :: MAX_ELEMENT_VERTICES =      6
  integer, parameter   :: MAX_SURFACE_VERTICES = 131072
  integer, parameter   :: MAX_SURFACE_ELEMENTS = 131072

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

    ! Is the mesh divided among processors (like in Front_Type)
    ! or is it shared among all of them (like in Surf_Type)
    logical :: mesh_divided

    ! Front to static Grid connectivity
    integer, allocatable :: cell_at_elem(:)
    integer, allocatable :: face_at_elem(:,:)

    contains

      procedure :: Allocate_Front
      procedure :: Calculate_Element_Centroids
      procedure :: Calculate_Element_Normals
      procedure :: Check_Elements
      procedure :: Clean_Front
      procedure :: Compress_Vertices
      procedure :: Compute_Distance_Function_And_Vof
      procedure :: Find_Front_Elements
      procedure :: Find_Sides
      procedure :: Find_Vertex_Elements
      procedure :: Handle_3_Points
      procedure :: Handle_4_Points
      procedure :: Handle_5_Points
      procedure :: Handle_6_Points
      procedure :: Initialize_Front
      procedure :: Mark_Cells_And_Faces
      procedure :: Place_Front_At_Value
      procedure :: Print_Front_Statistics
      ! procedure :: Save_Front

  end type

  contains

#   include "Front_Mod/Allocate_Front.f90"
#   include "Front_Mod/Calculate_Element_Centroids.f90"
#   include "Front_Mod/Calculate_Element_Normals.f90"
#   include "Front_Mod/Check_Elements.f90"
#   include "Front_Mod/Clean.f90"
#   include "Front_Mod/Compress_Vertices.f90"
#   include "Front_Mod/Compute_Distance_Function_And_Vof.f90"
#   include "Front_Mod/Find_Front_Elements.f90"
#   include "Front_Mod/Find_Sides.f90"
#   include "Front_Mod/Find_Vertex_Elements.f90"
#   include "Front_Mod/Handle_3_Points.f90"
#   include "Front_Mod/Handle_4_Points.f90"
#   include "Front_Mod/Handle_5_Points.f90"
#   include "Front_Mod/Handle_6_Points.f90"
#   include "Front_Mod/Initialize_Front.f90"
#   include "Front_Mod/Mark_Cells_And_Faces.f90"
#   include "Front_Mod/Place_Front_At_Value.f90"
#   include "Front_Mod/Print_Front_Statistics.f90"
!   include "Front_Mod/Save_Front.f90"

  end module
