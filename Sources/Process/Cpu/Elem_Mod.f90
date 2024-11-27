!==============================================================================!
  module Elem_Mod
!------------------------------------------------------------------------------!
!>  This module defines the Elem_Type, primarily in Front_Mod and Surf_Mod.
!>  Elem_Type includes attributes and functionalities tailored to the specific
!>  needs of these modules in relationship with meshing algorithms and even
!>  to sub-sequent use in Vof_Mod (curvatures, for example).  It is interesting
!>  to note that Elem_Mod is a direct descendant of a module with the same
!>  name which existed in EasyMesh, and the algorithms from EasyMesh have all
!>  been transferred to Surf_Mod, in almost one-to-one fashion.
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !---------------!
  !   Elem type   !
  !---------------!
  !> Encapsulates data (and an initialization function) needed for use of
  !> the Elem_Type objects in Front_Mod, Surf_Mod and even Vof_Mod, hence
  !> data which facilitates meshing algorithms and interface tracking.
  type Elem_Type

    integer :: nne         !! number of neighbouring element
    integer :: nv          !! number of vertices for the element
    integer :: ns          !! number of sides for the element
    integer :: v(24)       !! list of vertices
    integer :: e(24)       !! list of neighbouring elements
    integer :: s(24)       !! sides surrounding the element
    integer :: cell        !! cell at which the surface resides
    integer :: face        !! face which the surface intersects
    real    :: nx, ny, nz  !! surface normal vector
    real    :: xc, yc, zc  !! center of a sphere
    real    :: xe, ye, ze  !! center of element
    real    :: sx, sy, sz  !! surface area components
    real    :: area        !! element's surface area
    real    :: curv        !! curvature at the element

    contains
      procedure :: Initialize_Elem

  end type

  contains

#   include "Elem_Mod/Initialize_Elem.f90"

  end module
