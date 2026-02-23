!==============================================================================!
  module Face_Mod
!------------------------------------------------------------------------------!
!>  The Face_Mod module is designed for handling variables defined at the cell
!>  faces, complementing the cell-centered variables in T-Flows, which are
!>  handled by Var_Mod.
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!

  !---------------!
  !   Face type   !
  !---------------!
  !> The primary type in this module, Face_Type, is responsible for managing
  !> face-centered variables.  It is currently used to define volume fluxes at
  !> cell faces and it is unclear if its use will spread any further than that.
  type Face_Type

    type(Grid_Type), pointer :: pnt_grid  !! grid for which the
                                          !! variable is defined
    character(VL) :: name        !! variable name, always upper case and
                                 !! very short (4, defined in Const_Mod)
    real, allocatable :: n(:)    !! new value of the variable
    real, allocatable :: o(:)    !! old value of the variable
    real, allocatable :: oo(:)   !! older than old value
    real, allocatable :: avg(:)  !! average guessed value, guessed value

    ! Pointer to the device (Fortran view for OpenACC)
    real, pointer :: n_dev(:)

    ! CUDA device pointer handle
    type(c_ptr) :: n_dev_c

  end type

  contains

#   include "Face_Mod/Allocate.f90"

  end module
