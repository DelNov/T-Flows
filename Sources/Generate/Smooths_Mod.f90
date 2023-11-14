!==============================================================================!
  module Smooths_Mod
!------------------------------------------------------------------------------!
!>  The Smooths_Mod module in the Generate program provides data and
!>  functionality for smoothing a computational grid.  (I am not sure if
!>  smoothing is really involved from the main Generate's function.)
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !-----------------!
  !   Smooth_Type   !
  !-----------------!
  !> Encapsulates data for grid smoothing
  type Smooths_Type

    integer :: n_smooths    ! number of smoothing regions

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
