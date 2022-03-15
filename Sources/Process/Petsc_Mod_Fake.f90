!==============================================================================!
  module Petsc_Mod
!------------------------------------------------------------------------------!
!   This is a fake PETSc module, when code is compiled without it.             !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Native_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !----------------!
  !   Petsc type   !
  !----------------!
  type Petsc_Type

    type(Grid_Type), pointer :: pnt_grid

    ! Fake Petsc-related variables
    type(Matrix_Type) :: A         ! sparse matrix

    contains
      procedure :: Create_Petsc
      procedure :: Solve_Petsc

  end type

  logical, parameter :: PETSC_ACTIVE = .false.

  contains

    include 'Petsc_Mod_Fake/Create_Petsc.f90'
    include 'Petsc_Mod_Fake/Solve_Petsc.f90'

  end module
