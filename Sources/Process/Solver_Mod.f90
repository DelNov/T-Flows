#include "../Shared/Unused.h90"

!==============================================================================!
  module Solver_Mod
!------------------------------------------------------------------------------!
!>  The Solver_Mod module in T-Flows serves as a crucial abstraction layer for
!>  linear solvers within the software. It encapsulates the functionalities of
!>  both native solvers and PETSc solvers, offering a unified interface for
!>  their operations.  The design allows for easy toggling between different
!>  solver types without necessitating changes in the interaction with solvers
!>  in other parts of the code.  Despite the common purpose of solving linear
!>  systems, the implementation of native and PETSc solvers within this module
!>  is markedly different. This diversity renders Solver_Mod as one of the most
!>  heterogeneous modules in the entire T-Flows package.
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Var_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !-----------------!
  !   Solver type   !
  !-----------------!
  type Solver_Type

    type(Native_Type) :: Nat

    ! Linear solvers; NATIVE or PETSC
    integer :: solvers

    contains
      procedure :: Alias_Native
      procedure :: Create_Solver
      procedure :: End
      procedure :: Remove_Singular
      procedure :: Run
      procedure :: Set_Singular

  end type

  ! Linear solvers, native or PETSc
  integer, parameter :: NATIVE = 50021
  integer, parameter :: PETSC  = 50023

  contains

#   include "Solver_Mod/Alias_Native.f90"
#   include "Solver_Mod/Create_Solver.f90"
#   include "Solver_Mod/End.f90"
#   include "Solver_Mod/Linear_Solvers_Code.f90"
#   include "Solver_Mod/Remove_Singular.f90"
#   include "Solver_Mod/Run.f90"
#   include "Solver_Mod/Set_Singular.f90"

  end module
