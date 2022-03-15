#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>

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
  use Iso_C_Binding
  use PetscVec,   only: tVec
  use PetscMat,   only: tMat
  use PetscKSP,   only: tKSP, tPC
  use Native_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !----------------!
  !   Petsc type   !
  !----------------!
  type Petsc_Type

    type(Grid_Type), pointer :: pnt_grid

    ! Petsc-related variables
    type(tMat)                  :: A         ! sparse matrix
    type(tVec)                  :: x         ! solution vector
    type(tVec)                  :: b         ! right hand side
    type(tKSP)                  :: ksp       ! linear solver context
    type(tPc)                   :: pc        ! preconditioner
    type(PetscInt)              :: m_lower   ! unknowns in this proc.
    type(PetscInt)              :: m_upper   ! total number of unknowns
    type(PetscInt)              :: miter     ! maximum number of iterations
    type(PetscInt)              :: niter     ! performed number of iterations
    type(PetscInt), allocatable :: d_nnz(:)  ! diagonal stencil width per cell
    type(PetscInt), allocatable :: o_nnz(:)  ! off-diag. stencil width per cell
    type(PetscErrorCode)        :: err
    type(PetscEnum)             :: reason

    ! Global cell numbering for PETSc, which is
    ! different from T-Flows' and stars from zero   <---= IMPORTANT
    integer, allocatable :: glo(:)

    contains
      procedure :: Create_Petsc
      procedure :: Solve_Petsc

  end type

  logical, parameter :: PETSC_ACTIVE = .true.
  integer, parameter :: OUT_OF_ITS   = -3      ! KSP_DIVERGED_ITS

  contains

  include 'Petsc_Mod_True/Create_Petsc.f90'
  include 'Petsc_Mod_True/Solve_Petsc.f90'

  end module
