!==============================================================================!
  module Smooths_Mod
!------------------------------------------------------------------------------!
!   Type for smoothing a domain (a grid to be more correct).                   !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !-----------------!
  !   Smooth_Type   !
  !-----------------!
  type Smooths_Type

    integer              :: n_smooths    ! number of smoothing regions

    integer, allocatable :: iters(:)
    logical, allocatable :: in_x (:), in_y (:), in_z (:)
    real,    allocatable :: x_min(:), y_min(:), z_min(:)
    real,    allocatable :: x_max(:), y_max(:), z_max(:)
    real,    allocatable :: relax(:)

  end type

  !---------------------------!
  !   Member-like functions   !
  !---------------------------!
  contains

#   include "Smooths_Mod/Allocate_Smooths.f90"
#   include "Smooths_Mod/Grid.f90"

  end module
