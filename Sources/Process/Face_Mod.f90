!==============================================================================!
  module Face_Mod
!------------------------------------------------------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!

  !---------------!
  !   Face type   !
  !---------------!
  type Face_Type

    type(Grid_Type), pointer :: pnt_grid  ! grid for which it is defined

    character(VL) :: name        ! variable name, always upper case and
                                 ! very short (4, defined in Const_Mod)
    real, allocatable :: n(:)    ! new value
    real, allocatable :: o(:)    ! old value
    real, allocatable :: oo(:)   ! older than old value
    real, allocatable :: avg(:)  ! average guessed value, guessed value
  end type

  contains

#   include "Face_Mod/Allocate.f90"

  end module
