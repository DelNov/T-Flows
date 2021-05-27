!==============================================================================!
  module Surf_Mod
!------------------------------------------------------------------------------!
!   Module for Lagrangian particle tracking                                    !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Front_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !---------------!
  !   Surf type   !
  !---------------!
  type, extends(Front_Type) :: Surf_Type

    ! Logical array if cell has particles
    logical, allocatable :: cell_has_vertex(:)

    ! Working arrays, buffers for parallel version
    ! (Keyword "parameter: not allowed inside a type
    ! declaration. One might think of making a function)
    integer :: N_I_VARS = 3
    integer :: N_L_VARS = 2
    integer :: N_R_VARS = 8

    contains

      procedure :: Advance_Vertices
      procedure :: Allocate_Surf
      procedure :: Calculate_Curvatures_From_Edges
      procedure :: Calculate_Curvatures_From_Elems
      procedure :: Calculate_Curvatures_From_Verts
      procedure :: Find_Boundaries
      procedure :: Find_Surf_Elements
      procedure :: Improve_Mesh_Quality
      procedure :: Initialize_Surf
      procedure :: Place_Surf_At_Value
      procedure :: Refine
      procedure :: Relax_Geometry
      procedure :: Relax_Topology
      procedure :: Split_Element
      procedure :: Swap_Side
      procedure :: Smooth

  end type

  integer, allocatable :: i_work(:)
  logical, allocatable :: l_work(:)
  real,    allocatable :: r_work(:)

  contains

  include 'Surf_Mod/Advance_Vertices.f90'
  include 'Surf_Mod/Allocate_Surf.f90'
  include 'Surf_Mod/Clean.f90'
  include 'Surf_Mod/Calculate_Curvatures_From_Edges.f90'
  include 'Surf_Mod/Calculate_Curvatures_From_Elems.f90'
  include 'Surf_Mod/Calculate_Curvatures_From_Verts.f90'
  include 'Surf_Mod/Find_Boundaries.f90'
  include 'Surf_Mod/Find_Surf_Elements.f90'
  include 'Surf_Mod/Improve_Mesh_Quality.f90'
  include 'Surf_Mod/Initialize_Surf.f90'
  include 'Surf_Mod/Place_Surf_At_Value.f90'
  include 'Surf_Mod/Refine.f90'
  include 'Surf_Mod/Relax_Geometry.f90'
  include 'Surf_Mod/Relax_Topology.f90'
  include 'Surf_Mod/Split_Element.f90'
  include 'Surf_Mod/Swap_Side.f90'
  include 'Surf_Mod/Smooth.f90'
! include 'Surf_Mod/Bounce_Particle.f90'
! include 'Surf_Mod/Check_Periodicity.f90'
! include 'Surf_Mod/Exchange_Particles.f90'

  end module
