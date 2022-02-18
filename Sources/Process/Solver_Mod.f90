!==============================================================================!
  module Solver_Mod
!------------------------------------------------------------------------------!
!   This is a module which will entail native and PETSc solvers                !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Native_Mod
  use Petsc_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !-----------------!
  !   Solver type   !
  !-----------------!
  type Solver_Type

    type(Native_Type) :: Nat
    type(Petsc_Type)  :: Pet

    ! Linear solvers; native or PETSc
    integer :: solvers

    contains
      procedure :: Alias_Native
      procedure :: Create_Solver
      procedure :: Run

  end type

  ! Linear solvers, native or PETSc
  integer, parameter :: NATIVE = 50021
  integer, parameter :: PETSC  = 50023

  contains

  include 'Solver_Mod/Alias_Native.f90'
  include 'Solver_Mod/Create_Solver.f90'
  include 'Solver_Mod/Run.f90'

  end module
