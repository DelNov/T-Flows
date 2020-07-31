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

    character(VL) :: name      ! variable name, always upper case and
                               ! very short (4, defined in Const_Mod)
    real, allocatable :: n(:)           ! new value
    real, allocatable :: avg(:), star(:)! average guessed value, guessed value
    real, allocatable :: o(:), oo(:)    ! old and older then old
  end type

  contains

  include 'Face_Mod/Allocate_New_Only.f90'
  include 'Face_Mod/Allocate_New_And_Old.f90'

  end module
