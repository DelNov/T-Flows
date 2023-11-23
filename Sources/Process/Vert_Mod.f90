!==============================================================================!
  module Vert_Mod
!------------------------------------------------------------------------------!
!>  This module defines the Vert_Type, an extension of the Point_Type used
!>  primarily in Front_Mod and Surf_Mod. Vert_Type includes additional
!>  attributes and functionalities tailored to the specific needs of these
!>  modules, particularly in handling vertices on computational surfaces.
!------------------------------------------------------------------------------!
!   Features and Attributes                                                    !
!                                                                              !
!   * Extension of Point_Type: Inherits basic point properties and adds        !
!     specialized attributes for vertices.                                     !
!   * Coordinates history: Maintains both current and previous (old)           !
!     coordinates of the vertex.                                               !
!   * Neighbourhood information: Tracks neighbouring elements and vertices,    !
!     crucial for surface-related computations.                                !
!   * Boundary identification: Indicates whether the vertex is on a boundary.  !
!   * Curvature: Stores the curvature at the vertex, important for surface     !
!     processing.                                                              !
!   * State flags: Includes flags like 'deposited', 'escaped', and 'trapped'   !
!     to track the vertex's interaction with the computational domain.         !
!   * Smoothed values: Holds smoothed Volume of Fluid (VOF) values and their   !
!     gradients from the nearest cells, aiding in surface reconstruction.      !
!   * Initialization routine: Provides a method to initialize the vertex with  !
!     default or specified properties.                                         !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Point_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !---------------!
  !   Vert type   !
  !---------------!
  type, extends(Point_Type) :: Vert_Type
  !> Vert_Type is an extension of the Point_Type used in Front_Mod and
  !> Surf_Mod. Vert_Type includes additional attributes and functionalities
  !> such as handling vertices on computational surfaces.

    ! Old vertex's coordinates; new is in the parent
    real :: x_o  !! old x coordinate
    real :: y_o  !! old y coordinate
    real :: z_o  !! old z coordinate

    ! Needed in smoothing
    real :: sumx, sumy, sumz

    integer :: nne       !! number of neighbouring elements
    integer :: nnv       !! number of neighbouring vertices
    logical :: boundary  !! true if vertex is on a boundary
    integer :: e(24)     !! list of elements around the vertex

    ! The boundary face
    integer :: bnd_face  !! the closest boundary face

    ! Curvature at the vertex
    real :: curv  !! curvature at the vertex

    ! Needed for children, really
    logical :: deposited  !! deposited on a wall (particles)
    logical :: escaped    !! escaped from computational domain (particle)
    logical :: trapped    !! trapped on a surface (particle)

    ! Store values of the smoothed vof and
    ! its gradients from the nearest cells
    real :: smooth    !! store value of smoothed VOF
    real :: smooth_x  !! x gradient of the smoothed VOF
    real :: smooth_y  !! y gradient of the smoothed VOF
    real :: smooth_z  !! z gradient of the smoothed VOF

    ! Store closest cell coordinates
    real :: cell_x, cell_y, cell_z  !! closest cell coordinate

    contains
      procedure :: Initialize_Vert

  end type

  contains

#   include "Vert_Mod/Intialize_Vert.f90"

  end module
