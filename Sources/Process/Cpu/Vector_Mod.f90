!==============================================================================!
  module Vector_Mod
!------------------------------------------------------------------------------!
!>  The Vector_Mod module in T-Flows defines the structure and functionality
!>  for storage and handling of vector quantities associated with the
!>  computational grid.  The size of this vector is equal to the number of
!>  cells in the computational grid for which it is defined.
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!

  !-----------------!
  !   Vector type   !
  !-----------------!
  !> Designed for representing vector quantities defined over a grid
  !> and also provides a link to the computational grid for context.
  type Vector_Type

    type(Grid_Type), pointer :: pnt_grid  !! pointer to grid

    real, allocatable :: val(:)    !! values of the vector
  end type

  contains

#   include "Vector_Mod/Allocate.f90"

  end module
