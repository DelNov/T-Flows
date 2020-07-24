!==============================================================================!
  module Surf_Mod
!------------------------------------------------------------------------------!
!   Module for Lagrangian particle tracking                                    !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Math_Mod
  use Sort_Mod
  use Comm_Mod
  use Grid_Mod
  use Var_Mod
  use Field_Mod
  use Solver_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !--------------------------------!
  !   A few important parameters   !
  !--------------------------------!
  integer, parameter   :: MAX_ELEMENT_VERTICES =      6
  integer, parameter   :: MAX_SURFACE_VERTICES = 131072
  integer, parameter   :: MAX_SURFACE_ELEMENTS = 131072

  !---------------!
  !   Vert type   !
  !---------------!
  type Vert_Type

    ! Vertex's coordinates; new and old
    real :: x_n, x_o
    real :: y_n, y_o
    real :: z_n, z_o

    ! Vertex's velocity (if needed)
    real :: u
    real :: v
    real :: w

    real :: sumx, sumy, sumz

    integer :: nne       ! number of neighbouring elements
    integer :: nnv       ! number of neighbouring vertices
    logical :: boundary  ! is vertex on a boundary

    ! The closest cell, node, boundary cell and face
    integer :: cell
    integer :: node
    integer :: bnd_cell
    integer :: bnd_face

    ! Vertex departure from domain 
    logical :: escaped

    ! Vertex inside the subdomain
    integer :: proc
    integer :: buff

    real :: curv
  end type

  !---------------!
  !   Elem type   !
  !---------------!
  type Elem_Type

    integer :: nne         ! number of neighbouring element
    integer ::  i,  j,  k
    integer :: ei, ej, ek
    integer :: si, sj, sk
    real    :: nx, ny, nz  ! surface normal vector
    real    :: area
    real    :: curv

  end type

  !---------------!
  !   Side type   !
  !---------------!
  type Side_Type

    integer :: ei, ea, eb      ! element undefined, elements left and right
    integer :: a, b, c, d
    real    :: length
    logical :: boundary

  end type

  !---------------!
  !   Surf type   !
  !---------------!
  type Surf_Type

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
  end type

  ! Working arrays, buffers for parallel version
  integer, parameter   :: N_I_VARS = 3
  integer, parameter   :: N_L_VARS = 2
  integer, parameter   :: N_R_VARS = 8
  integer, allocatable :: i_work(:)
  logical, allocatable :: l_work(:)
  real,    allocatable :: r_work(:)

  contains

! include 'Surf_Mod/Advance_Particles.f90'
  include 'Surf_Mod/Allocate.f90'
  include 'Surf_Mod/Calculate_Element_Normals.f90'
  include 'Surf_Mod/Clean.f90'
  include 'Surf_Mod/Count_Elements_Neighbours.f90'
  include 'Surf_Mod/Count_Vertex_Elements.f90'
  include 'Surf_Mod/Compress_Vertices.f90'
  include 'Surf_Mod/Calculate_Curvatures_From_Edges.f90'
  include 'Surf_Mod/Calculate_Curvatures_From_Elems.f90'
  include 'Surf_Mod/Calculate_Curvatures_From_Verts.f90'
  include 'Surf_Mod/Find_Boundaries.f90'
  include 'Surf_Mod/Find_Sides.f90'
  include 'Surf_Mod/Find_Nearest_Cell.f90'
  include 'Surf_Mod/Find_Nearest_Node.f90'
  include 'Surf_Mod/Handle_3_Points.f90'
  include 'Surf_Mod/Handle_4_Points.f90'
  include 'Surf_Mod/Handle_5_Points.f90'
  include 'Surf_Mod/Handle_6_Points.f90'
  include 'Surf_Mod/Place_At_Var_Value.f90'
  include 'Surf_Mod/Print_Statistics.f90'
  include 'Surf_Mod/Refine.f90'
  include 'Surf_Mod/Relax_Geometry.f90'
  include 'Surf_Mod/Relax_Topology.f90'
  include 'Surf_Mod/Split_Element.f90'
  include 'Surf_Mod/Swap_Side.f90'
  include 'Surf_Mod/Smooth.f90'
! include 'Surf_Mod/Bounce_Particle.f90'
! include 'Surf_Mod/Check_Periodicity.f90'
! include 'Surf_Mod/Exchange_Particles.f90'
! include 'Surf_Mod/Move_Particle.f90'

  end module
