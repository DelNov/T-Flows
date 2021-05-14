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

    type(Grid_Type),  pointer :: pnt_grid  ! grid for which it is defined
    type(Field_Type), pointer :: pnt_flow  ! flow field for which it is defined

    integer                      :: n_elems
    integer                      :: n_verts
    integer                      :: n_sides
    type(Vert_Type), allocatable :: vert(:)
    type(Elem_Type), allocatable :: elem(:)
    type(Side_Type), allocatable :: side(:)

    ! Front to static grid connectivity
    integer, allocatable :: cell_at_elem(:)
    integer, allocatable :: face_at_elem(:,:)

    contains

      procedure :: Allocate_Front
      procedure :: Calculate_Element_Centroids
      procedure :: Calculate_Element_Normals
      procedure :: Clean_Front
      procedure :: Compress_Vertices
      procedure :: Compute_Distance_Function_And_Vof
      procedure :: Calculate_Curvatures_From_Elems
      procedure :: Find_Connectivity
      procedure :: Find_Nearest_Cell
      procedure :: Find_Nearest_Node
      procedure :: Find_Vertex_Elements
      procedure :: Handle_3_Points
      procedure :: Handle_4_Points
      procedure :: Handle_5_Points
      procedure :: Handle_6_Points
      procedure :: Initialize_Front
      procedure :: Mark_Cells_And_Faces
      procedure :: Place_At_Var_Value
      procedure :: Print_Statistics
      procedure :: Save_Front

  end type

  contains

  include 'Front_Mod/Allocate.f90'
  include 'Front_Mod/Calculate_Element_Centroids.f90'
  include 'Front_Mod/Calculate_Element_Normals.f90'
  include 'Front_Mod/Clean.f90'
  include 'Front_Mod/Compress_Vertices.f90'
  include 'Front_Mod/Compute_Distance_Function_And_Vof.f90'
  include 'Front_Mod/Calculate_Curvatures_From_Elems.f90'
  include 'Front_Mod/Find_Connectivity.f90'
  include 'Front_Mod/Find_Nearest_Cell.f90'
  include 'Front_Mod/Find_Nearest_Node.f90'
  include 'Front_Mod/Find_Vertex_Elements.f90'
  include 'Front_Mod/Handle_3_Points.f90'
  include 'Front_Mod/Handle_4_Points.f90'
  include 'Front_Mod/Handle_5_Points.f90'
  include 'Front_Mod/Handle_6_Points.f90'
  include 'Front_Mod/Initialize.f90'
  include 'Front_Mod/Mark_Cells_And_Faces.f90'
  include 'Front_Mod/Place_At_Var_Value.f90'
  include 'Front_Mod/Print_Statistics.f90'
  include 'Front_Mod/Save.f90'

  end module
