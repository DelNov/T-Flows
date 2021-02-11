!==============================================================================!
  module Elem_Mod
!------------------------------------------------------------------------------!
!   Storage for Elem_Type used by Front_Mod and Surf_mod                       !
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !---------------!
  !   Elem type   !
  !---------------!
  type Elem_Type

    integer :: nne         ! number of neighbouring element
    integer ::  i,  j,  k
    integer :: ei, ej, ek
    integer :: si, sj, sk
    real    :: nx, ny, nz  ! surface normal vector
    real    :: xc, yc, zc  ! center of a sphere
    real    :: xe, ye, ze  ! center of element
    real    :: area
    real    :: curv

  end type

  end module
