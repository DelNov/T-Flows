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
    integer, allocatable :: buff_i(:), buff_j(:), buff_k(:), buff_n(:)
    real,    allocatable :: buff_x(:), buff_y(:), buff_z(:), buff_v(:)

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
      procedure, private :: Distribute_Mesh
      procedure, private :: Distribute_Smooth
      procedure, private :: Distribute_Cell_Coords
      procedure, private :: Refine
      procedure, private :: Relax_Geometry
      procedure, private :: Relax_Topology
      procedure, private :: Smooth_Surf
      procedure, private :: Split_Element
      procedure, private :: Swap_Side

  end type

  contains

#   include "Surf_Mod/Advance_Vertices.f90"
#   include "Surf_Mod/Allocate_Surf.f90"
#   include "Surf_Mod/Clean.f90"
#   include "Surf_Mod/Calculate_Curvatures_From_Edges.f90"
#   include "Surf_Mod/Calculate_Curvatures_From_Elems.f90"
#   include "Surf_Mod/Calculate_Curvatures_From_Verts.f90"
#   include "Surf_Mod/Find_Boundaries.f90"
#   include "Surf_Mod/Find_Surf_Elements.f90"
#   include "Surf_Mod/Improve_Mesh_Quality.f90"
#   include "Surf_Mod/Initialize_Surf.f90"
#   include "Surf_Mod/Place_Surf_At_Value.f90"
#   include "Surf_Mod/Distribute_Mesh.f90"
#   include "Surf_Mod/Distribute_Smooth.f90"
#   include "Surf_Mod/Distribute_Cell_Coords.f90"
#   include "Surf_Mod/Refine.f90"
#   include "Surf_Mod/Relax_Geometry.f90"
#   include "Surf_Mod/Relax_Topology.f90"
#   include "Surf_Mod/Smooth_Surf.f90"
#   include "Surf_Mod/Split_Element.f90"
#   include "Surf_Mod/Swap_Side.f90"
!   include "Surf_Mod/Bounce_Particle.f90"
!   include "Surf_Mod/Check_Periodicity.f90"

  end module
