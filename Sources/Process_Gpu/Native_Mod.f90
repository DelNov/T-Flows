#include "../Shared/Browse.h90"

!==============================================================================!
  module Native_Mod
!------------------------------------------------------------------------------!
  use Work_Mod
  use Linalg_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  integer, parameter :: MATRIX_UVW  =  1
  integer, parameter :: MATRIX_PP   =  2
  integer, parameter :: MATRIX_T    =  3
  integer, parameter :: MATRIX_ONE  =  1

  !-----------------!
  !   Native type   !
  !-----------------!
  type Native_Type

    type(Grid_Type), pointer :: pnt_grid  !! pointer to the numerical grid

    ! Connectivity matrix
    type(Sparse_Con_Type) :: C  !! connectivity matrix for all variables

    ! Matrices for all variables.  In this early stage of development,
    ! index 1 will belong to momentum, index 2 to pressure
    type(Sparse_Val_Type) :: A(3)  !! system value matrices

    ! If the following data member is set to .true., all variables will
    ! be discretized in one memory space, one matrix.  Economic memory-
    ! wise, but more time is spent on forming matrices.  If, however,
    ! it is set to .false., each variable will have its own memory space
    ! for its matrix.  That will use more memory, but will be faster.
    logical :: use_one_matrix  !! matrix space reused?

    ! Right-hand side vector for all variables
    real, allocatable :: b(:)

    ! Vectors used with native solvers
    real, allocatable :: p(:)      !! helping vector
    real, allocatable :: q(:)      !! helping vector
    real, allocatable :: r(:)      !! residual vector

    contains
      procedure :: Cg             !! conjugate gradient solver
      procedure :: Create_Native  !! creates native solver context

  end type

  contains

#   include "Native_Mod/Cg.f90"
#   include "Native_Mod/Create_Native.f90"

  end module
