!==============================================================================!
  module Vert_Mod
!------------------------------------------------------------------------------!
!   Storage for Vert_Type used by Front_Mod and Surf_mod                       !
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !---------------!
  !   Vert type   !
  !---------------!
  type Vert_Type

    ! Vertex's coordinates; new and old
    real :: x_n, x_o
    real :: y_n, y_o
    real :: z_n, z_o

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
    integer :: cell
    integer :: node
    integer :: bnd_cell
    integer :: bnd_face

    ! Vertex departure from domain 
    logical :: escaped

    ! Vertex inside the subdomain
    integer :: proc
    integer :: buff

    ! Curvature at the vertex
    real :: curv
  end type

  end module
