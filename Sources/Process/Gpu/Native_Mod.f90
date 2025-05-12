#include "../../Shared/Browse.h90"
#include "../../Shared/Unused.h90"

!==============================================================================!
  module Native_Mod
!------------------------------------------------------------------------------!
  use Work_Mod
  use Linalg_Mod
  use Gpu_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !-----------------!
  !   Native type   !
  !-----------------!
  type Native_Type

    type(Grid_Type), pointer :: pnt_grid  !! pointer to the numerical grid

    ! Connectivity matrix
    type(Sparse_Con_Type) :: C  !! connectivity matrix for all variables

    ! Matrix for all variables.  (It used to be array, but I abandoned it)
    type(Sparse_Val_Type) :: A  !! system value matrices

    ! Right-hand side vector for all variables
    real, allocatable :: b(:)

    ! Vectors used with native solvers
    real, allocatable :: p(:)      !! helping vector
    real, allocatable :: q(:)      !! helping vector
    real, allocatable :: r(:)      !! residual vector

    contains
      procedure :: Cg             !! conjugate gradient solver
      procedure :: Create_Native  !! creates native solver context

      ! Procedures to copy native solver to device
      procedure :: Copy_Native_To_Device
      procedure :: Destroy_Native_On_Device

  end type

  contains

#   include "Native_Mod/Cg.f90"
#   include "Native_Mod/Create_Native.f90"

    ! Procedures to copy native solver to device
#   include "Native_Mod/Gpu/Copy_To_Device.f90"
#   include "Native_Mod/Gpu/Destroy_On_Device.f90"


  end module
