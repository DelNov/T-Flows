!==============================================================================!
  subroutine Initialize_Vert(Vert, Grid)
!------------------------------------------------------------------------------!
!>  This subroutine is responsible for initializing the Vert_Type instance.
!>  It sets up all the necessary attributes for a vertex in the computational
!>  grid, providing a base state from which the vertex's properties can evolve
!>  during the simulation process.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Inheritance initialization: Calls the constructor of the parent          !
!     Point_Type to initialize basic point properties.                         !
!   * Coordinate setup: Initializes both current and previous coordinates      !
!     of the vertex.                                                           !
!   * Smoothing components: Sets initial values for smoothing-related          !
!     variables, likely used in surface reconstruction or processing.          !
!   * Neighborhood counters: Initializes the counters for neighboring          !
!     elements and vertices, crucial for interactions and computations         !
!     involving the vertex.                                                    !
!   * Boundary flag: Establishes the initial boundary state of the vertex.     !
!   * Element list: Initializes the list of elements surrounding the vertex,   !
!     important for mesh-based calculations.                                   !
!   * Closest features: Sets initial values for the nearest cell, node,        !
!     boundary cell, and face.                                                 !
!   * Curvature: Initializes the curvature at the vertex, a key geometric      !
!     property.                                                                !
!   * State flags: Sets initial states for 'deposited', 'escaped', and         !
!     'trapped' flags, reflecting the vertex's interaction with the domain.    !
!   * Smoothed values: Prepares variables for storing smoothed values and      !
!     their gradients, aiding in calculations involving the vertex.            !
!   * Cell coordinates: Initializes coordinates for the closest cell,          !
!     providing a reference for spatial calculations.                          !
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vert_Type)        :: Vert  !! vertex being initialized
  type(Grid_Type), target :: Grid  !! grid on which the vertex is initialized
!==============================================================================!

  ! Call parent's constructor
  call Vert % Initialize_Point(Grid)

  ! Old vertex's coordinates; new is in the parent
  Vert % x_o = 0.0
  Vert % y_o = 0.0
  Vert % z_o = 0.0

  ! I forgot what these were, something for smoothing, I think
  Vert % sumx = 0.0
  Vert % sumy = 0.0
  Vert % sumz = 0.0

  Vert % nne      = 0        ! number of neighbouring elements
  Vert % nnv      = 0        ! number of neighbouring vertices
  Vert % boundary = .false.  ! is vertex on a boundary
  Vert % e(1:24)  = 0        ! list of elements around the vertex

  ! The closest cell, node, boundary cell and face
  Vert % bnd_face = 0

  ! Curvature at the vertex
  Vert % curv = 0.0

  Vert % deposited = .false. ! deposited on a wall (particles)
  Vert % escaped   = .false. ! escaped from computational domain
  Vert % trapped   = .false. ! trapped on a surface

  Vert % smooth   = 0.0
  Vert % smooth_x = 0.0
  Vert % smooth_y = 0.0
  Vert % smooth_z = 0.0

  Vert % cell_x = 0.0
  Vert % cell_y = 0.0
  Vert % cell_z = 0.0

  end subroutine
