!==============================================================================!
  module Point_Mod
!------------------------------------------------------------------------------!
!>  This module is responsible for representing the concept of a point, the
!>  most essential geometrical element in the process of grid generation.
!>  Together with Line_Mod (Line_Type) and Block_Mod (Block_Type) it helps to
!>  define a computational domain in Generate's Domain_Mod.
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !----------------!
  !   Point type   !
  !----------------!
  !> Encapsulates the properties of a point,
  !> which boil down to three coordinates.
  type Point_Type

    ! Point coordinates
    real :: x  !! x cooridnate
    real :: y  !! y coordinate
    real :: z  !! z coordinate

  end type

  end module
