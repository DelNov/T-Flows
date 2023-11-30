#include "../Shared/Assert.h90"
#include "../Shared/Browse.h90"
#include "../Shared/Unused.h90"

!==============================================================================!
  module Surf_Mod
!------------------------------------------------------------------------------!
!>  The Surf_Mod module in T-Flows is an advanced alternative to Front_Mod
!>  (from which it derived), designed for generating high-quality mesh surfaces
!>  independent of the underlying computational grid.  This independence makes
!>  Surf_Mod more versatile and capable of creating more refined meshes,
!>  especially beneficial in simulations involving complex interfaces.
!>  However, it has its share of shortcomings, such as more challenging
!>  parallelization and incompatibility with background polyhedral grids.
!>  Surf_Mod, derived from Front_Type, includes functionalities such as mesh
!>  refinement, smoothing, and handling complex geometrical points. These
!>  enable precise representation of surfaces defined by scalar fields (VOF).
!>  It is worth noting that Surf_Mod works with triangular grids only, and
!>  that the algorithms for improving the mesh quality all come from TRIPOS
!>  (https://github.com/Niceno/TRIPOS)
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Front_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !---------------!
  !   Surf type   !
  !---------------!
  !> he Surf_Type, an extension of Front_Type, focuses on the manipulation of
  !> triangular grids, facilitating high-quality mesh generation independent
  !> of the underlying computational grid. This type includes specialized
  !> functionalities for mesh improvement, mesh smoothing and curvatures
  !> calculation.
  type, extends(Front_Type) :: Surf_Type

    ! Logical array if cell has particles
    logical, allocatable :: cell_has_vertex(:)

    ! Working arrays, buffers for parallel version
    integer, allocatable :: buff_i(:), buff_j(:), buff_k(:), buff_n(:)
    real,    allocatable :: buff_x(:), buff_y(:), buff_z(:), buff_v(:)

    contains

      procedure          :: Advance_Vertices
      procedure          :: Allocate_Surf
      procedure          :: Calculate_Curvatures_From_Edges
      procedure          :: Calculate_Curvatures_From_Elems
      procedure          :: Calculate_Curvatures_From_Verts
      procedure          :: Compress_Surf_Vertices
      procedure          :: Find_Boundaries
      procedure          :: Find_Surf_Elements
      procedure          :: Handle_3_Points
      procedure          :: Handle_4_Points
      procedure          :: Handle_5_Points
      procedure          :: Handle_6_Points
      procedure          :: Improve_Mesh_Quality
      procedure          :: Initialize_Surf
      procedure          :: Place_Surf_At_Value
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
#   include "Surf_Mod/Compress_Surf_Vertices.f90"
#   include "Surf_Mod/Find_Boundaries.f90"
#   include "Surf_Mod/Find_Surf_Elements.f90"
#   include "Surf_Mod/Handle_3_Points.f90"
#   include "Surf_Mod/Handle_4_Points.f90"
#   include "Surf_Mod/Handle_5_Points.f90"
#   include "Surf_Mod/Handle_6_Points.f90"
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
