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

    ! Logical array if cell has particles
    logical, allocatable :: cell_has_vertex(:)

    ! Working arrays, buffers for parallel version
    ! (Keyword "parameter: not allowed inside a type
    ! declaration. One might think of making a function)
    integer :: N_I_VARS = 3
    integer :: N_L_VARS = 2
    integer :: N_R_VARS = 8
  end type

  integer, allocatable :: i_work(:)
  logical, allocatable :: l_work(:)
  real,    allocatable :: r_work(:)

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
  include 'Front_Mod/Place_At_Var_Value.f90'
  include 'Front_Mod/Print_Statistics.f90'

  end module
