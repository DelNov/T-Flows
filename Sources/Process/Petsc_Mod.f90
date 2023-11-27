#include "../Shared/Assert.h90"
#include "../Shared/Browse.h90"
#include "../Shared/Unused.h90"

#if T_FLOWS_PETSC == 1
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>
#endif

!==============================================================================!
  module Petsc_Mod
!------------------------------------------------------------------------------!
!>  The Petsc_Mod module in T-Flows is an essential component for integrating
!>  PETSc (Portable, Extensible Toolkit for Scientific Computation) solvers
!>  into the simulation environment.  This module is designed to handle linear
!>  algebra operations, particularly those relevant to solving sparse linear
!>  systems, arising after discretization of partial differential equations
!>  with T-Flows.
!------------------------------------------------------------------------------!
!   Conditional Compilation                                                    !
!                                                                              !
!   * The module's functionality is contingent on the T_FLOWS_PETSC flag.      !
!     When PETSc is enabled, the module includes PETSc-specific data types     !
!     and operations; otherwise, it resorts to placeholder or 'fake'           !
!     implementations.
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
# if T_FLOWS_PETSC == 1
    use Iso_C_Binding
    use PetscVec,   only: tVec
    use PetscMat,   only: tMat
    use PetscKSP,   only: tKSP, tPC
# endif
  use Native_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  ! Maximum number of PETSc objects to create
  integer, parameter :: MAX_PETSC_MEMBERS = 32  !! maximum number of PETSc
                                                !! instances to create
  !----------------!
  !   Petsc type   !
  !----------------!
  !> Encapsulates data and procedures which manage the interface between
  !> T-Flows and PETSc, enabling the use of PETSc's robust linear solvers
  !> and preconditioners for solving sparse matrix equations.
  type Petsc_Type

    type(Grid_Type), pointer :: pnt_grid

#   if T_FLOWS_PETSC == 1
      ! Petsc-related variables
      type(tMat)           :: A         !! sparse matrix in PETSc format
      type(tVec)           :: x         !! solution vector
      type(tVec)           :: b         !! right hand side (source) vector
      type(tKSP)           :: ksp       !! linear solver context
      type(tPc)            :: pc        !! preconditioner
      type(PetscInt)       :: m_lower   !! number of unknowns in this processor
      type(PetscInt)       :: m_upper   !! total number of unknowns
      type(PetscInt)       :: miter     !! maximum number of iterations
      type(PetscInt)       :: niter     !! performed number of iterations
      type(PetscErrorCode) :: err       !! PETSc error code
      type(PetscEnum)      :: reason    !! reason for exiting a linear sovler
      logical              :: matrix_coppied = .false.  !! parameter which is
        !! not a part of PETSc library, but is introduced to mitigiate
        !! excessive copying of matrices from T-Flow space to PETSc space
      logical              :: precond_formed = .false.  !! paramater which is
        !! not a part of PETSc library, but is introduced to mitigiate
        !! too frequent formation of preconditioning matrices in PETSc.
#   else
      ! Fake Petsc-related variables
      type(Matrix_Type) :: A  !! sparse matrix as defined in T-Flows
#   endif

    contains
      procedure :: Create_Petsc
      procedure :: Destroy_Petsc
      procedure :: Solve_Petsc

  end type

# if T_FLOWS_PETSC == 1
    logical, parameter :: PETSC_ACTIVE = .true.   !! is PETSc active?
    integer, parameter :: OUT_OF_ITS   = -3       !! alias for KSP_DIVERGED_ITS
# else
    logical, parameter :: PETSC_ACTIVE = .false.  !! is PETSc active?
# endif

    !---------------------------------------------!
    !   I couldn't find a better way to do this   !
    !---------------------------------------------!
    logical       :: petsc_is_reporting = .false.     !! is PETSc reporting?
    character(SL) :: petsc_options(MAX_STRING_ITEMS)  !! options for PETSc

    !---------------------!
    !   Work PETSC type   !
    !---------------------!
    !> A pool of PETSc instances used in a particular simulation
    type Work_Petsc_Type
      integer          :: n_members = 0              !! number of instances
      type(Petsc_Type) :: Member(MAX_PETSC_MEMBERS)  !! list of members
    end type

    type(Work_Petsc_Type), target :: Work_Pet  !! object of Work_Petsc_Type

  contains

# if T_FLOWS_PETSC == 1
#   include "Petsc_Mod/True/Create_Petsc.f90"
#   include "Petsc_Mod/True/Destroy_Petsc.f90"
#   include "Petsc_Mod/True/Solve_Petsc.f90"
# else
#   include "Petsc_Mod/Fake/Create_Petsc.f90"
#   include "Petsc_Mod/Fake/Destroy_Petsc.f90"
#   include "Petsc_Mod/Fake/Solve_Petsc.f90"
# endif

  end module
