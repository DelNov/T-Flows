!==============================================================================!
  module Vert_Mod
!------------------------------------------------------------------------------!
!   Storage for Vert_Type used by Front_Mod and Surf_mod                       !
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

    ! Old vertex's coordinates; new is in the parent
    real :: x_o
    real :: y_o
    real :: z_o

    ! Needed in smoothing
    real :: sumx, sumy, sumz

    integer :: nne       ! number of neighbouring elements
    integer :: nnv       ! number of neighbouring vertices
    logical :: boundary  ! is vertex on a boundary
    integer :: e(24)     ! list of elements around the vertex

    ! The closest cell, node, boundary cell and face
    integer :: bnd_face

    ! Curvature at the vertex
    real :: curv

    ! Needed for children, really
    logical :: deposited  ! deposited on a wall (particles)
    logical :: escaped    ! escaped from computational domain
    logical :: trapped    ! trapped on a surface

    ! Store values of the smoothed vof and
    ! its gradients from the nearest cells
    real :: smooth, smooth_x, smooth_y, smooth_z

    ! Store closest cell coordinates
    real :: cell_x, cell_y, cell_z

    contains
      procedure :: Initialize_Vert

  end type

  contains

  include 'Vert_Mod/Intialize_Vert.f90'

  end module
