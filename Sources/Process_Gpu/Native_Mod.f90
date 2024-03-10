!==============================================================================!
  module Native_Mod
!------------------------------------------------------------------------------!
  use Linalg_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !-----------------!
  !   Native type   !
  !-----------------!
  type Native_Type

    type(Grid_Type), pointer :: pnt_grid  !! pointer to the numerical grid

    ! Matrix for all variables except momentum
    type(Sparse_Type) :: A  !! system matrix for all variables except momentum

    ! Matrix for momentum equations
    type(Sparse_Type) :: M  !! system matrix for all variables except momentum

    ! Right-hand side vector for all variables
    real, allocatable :: b(:)

    ! Vectors used with native solvers
    real, allocatable :: p(:)      !! helping vector
    real, allocatable :: q(:)      !! helping vector
    real, allocatable :: r(:)      !! residual vector

    contains
      procedure :: Cg             !! conjugate gradient solver
      procedure :: Create_Native  !! creates native solver context
      procedure :: Prec_Form      !! forms a preconditioning matrix

  end type

  contains

#   include "Native_Mod/Cg.f90"
#   include "Native_Mod/Create_Native.f90"
#   include "Native_Mod/Prec_Form.f90"

  end module
