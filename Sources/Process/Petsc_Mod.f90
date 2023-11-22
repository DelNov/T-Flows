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
!   Module used for PETSc linear solvers.                                      !
!                                                                              !
!   Note: This module has all member procedures bundled in one file.  When I   !
!         included them from separate files, I experienced difficulties with   !
!         include files above and definition of PETSc types.  Go figure?       !
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
  integer, parameter :: MAX_PETSC_MEMBERS = 32

  !----------------!
  !   Petsc type   !
  !----------------!
  type Petsc_Type

    type(Grid_Type), pointer :: pnt_grid

#   if T_FLOWS_PETSC == 1
      ! Petsc-related variables
      type(tMat)           :: A         ! sparse matrix
      type(tVec)           :: x         ! solution vector
      type(tVec)           :: b         ! right hand side
      type(tKSP)           :: ksp       ! linear solver context
      type(tPc)            :: pc        ! preconditioner
      type(PetscInt)       :: m_lower   ! unknowns in this proc.
      type(PetscInt)       :: m_upper   ! total number of unknowns
      type(PetscInt)       :: miter     ! maximum number of iterations
      type(PetscInt)       :: niter     ! performed number of iterations
      type(PetscErrorCode) :: err
      type(PetscEnum)      :: reason
      logical              :: matrix_coppied = .false.
      logical              :: precond_formed = .false.
#   else
      ! Fake Petsc-related variables
      type(Matrix_Type) :: A         ! sparse matrix
#   endif

    contains
      procedure :: Create_Petsc
      procedure :: Destroy_Petsc
      procedure :: Solve_Petsc

  end type

# if T_FLOWS_PETSC == 1
    logical, parameter :: PETSC_ACTIVE = .true.
    integer, parameter :: OUT_OF_ITS   = -3      ! KSP_DIVERGED_ITS
# else
    logical, parameter :: PETSC_ACTIVE = .false.
# endif

    !---------------------------------------------!
    !   I couldn't find a better way to do this   !
    !---------------------------------------------!
    logical       :: petsc_is_reporting = .false.     ! is PETSc reporting
    character(SL) :: petsc_options(MAX_STRING_ITEMS)  ! options for PETSc

    !---------------------!
    !   Work PETSC type   !
    !---------------------!
    type Work_Petsc_Type
      integer          :: n_members = 0
      type(Petsc_Type) :: Member(MAX_PETSC_MEMBERS)
    end type

    type(Work_Petsc_Type), target :: Work_Pet

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
