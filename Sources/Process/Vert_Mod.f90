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

    ! Vertex's new velocity (if ever needed)
    real :: u_n
    real :: v_n
    real :: w_n

    real :: sumx, sumy, sumz

    integer :: nne       ! number of neighbouring elements
    integer :: nnv       ! number of neighbouring vertices
    logical :: boundary  ! is vertex on a boundary
    integer :: e(24)     ! list of elements around the vertex

    ! The closest cell, node, boundary cell and face
    integer :: bnd_face

    ! Vertex departure from domain
    logical :: escaped

    ! Curvature at the vertex
    real :: curv

    contains
      procedure :: Initialize_Vert

  end type

  contains

  include 'Vert_Mod/Intialize_Vert.f90'

  end module
