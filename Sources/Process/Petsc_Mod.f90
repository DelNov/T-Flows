!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!                                                                              !
!   PETSc is installed on the system -> take true PETSc module                 !
!                                                                              !
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
#ifdef PETSC_DIR

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
!         include files above.                                                 !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Iso_C_Binding
  use Solver_Mod
  use PetscVec,   only: tVec
  use PetscMat,   only: tMat
  use PetscKSP,   only: tKSP, tPC
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  type(PetscInt) :: PETSC_ONE = 1  ! one

  !------------------!
  !   Solvers type   !
  !------------------!
  type Petsc_Type

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
    type(PetscInt), allocatable :: col(:)    ! column indices
    type(PetscErrorCode)        :: err

    contains
      procedure :: Create_Petsc
      procedure :: Solve

  end type

  contains

!==============================================================================!
  subroutine Create_Petsc(Pet, Sol, Grid)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Petsc_Type)        :: Pet
  type(Solver_Type)        :: Sol
  type(Grid_Type),  target :: Grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, c1, c2, s, j, nonzeros
!==============================================================================!

  if(this_proc < 2) print *, '# Determining matrix topology.'

  !----------------------!
  !   Initialize PETSc   !
  !----------------------!
  call C_Petsc_Initialize()

  !-------------------------!
  !   Set up PETSc matrix   !
  !-------------------------!

  ! Total number of unknowns and unknowns in this processor only
  Pet % m_lower = Grid % n_cells - Grid % comm % n_buff_cells
  Pet % m_upper = Grid % comm % nc_tot

  ! Allocate memory for array with number of non-zero entries per row
  allocate(Pet % d_nnz(Grid % n_cells))
  Pet % d_nnz(:) = 1

  ! Compute stencil widths
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)
    if(c2 > 0) then
      Pet % d_nnz(c1) = Pet % d_nnz(c1) + 1
      Pet % d_nnz(c2) = Pet % d_nnz(c2) + 1
    end if
  end do

  ! Create PETSc matrix
  call C_Petsc_Mat_Create_Seq_Aij(Pet % A, Pet % m_lower, Pet % d_nnz)

  !-----------------------------------------!
  !   Set column indices for PETSc matrix   !
  !-----------------------------------------!

  ! Work out number of non-zero entries
  nonzeros = sum(Pet % d_nnz(1:Pet % m_lower))

  ! Allocate memory for it
  allocate(Pet % col(nonzeros));  Pet % col(:) = 0

  ! For the time being, copy the values from T-Flows minus 1
  Pet % col(1:nonzeros) = Sol % A % col(1:nonzeros) - PETSC_ONE

  ! Set column indices (now this is really tough)
  call C_Petsc_Mat_Seq_Aij_Set_Column_Indices(Pet % A, Pet % col)

  !--------------------------!
  !   Create PETSc vectors   !
  !--------------------------!
  call C_Petsc_Vec_Create_Seq(Pet % x, Pet % m_lower)
  call C_Petsc_Vec_Create_Seq(Pet % b, Pet % m_lower)

  !-------------------------!
  !   Create PETSc solver   !
  !-------------------------!
  call C_Petsc_Ksp_Create(Pet % ksp)

  if(this_proc < 2) print *, '# Finished !'

  end subroutine

!==============================================================================!
  subroutine Solve(Pet, Sol, A, x, b, miter, niter, tol)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Petsc_Type)    :: Pet
  type(Solver_Type)    :: Sol
  type(Matrix_Type)    :: A
  real                 :: x(-Sol % pnt_grid % n_bnd_cells :  &
                             Sol % pnt_grid % n_cells)
  real                 :: b( Sol % pnt_grid % n_cells)
  integer, intent(in)  :: miter
  integer, intent(out) :: niter
  real,    intent(in)  :: tol
!-----------------------------------[Locals]-----------------------------------!
  type(PetscInt) :: i, j, k
  character(SL)  :: sols           ! fortran string to store solver
  character(SL)  :: precs          ! fortran string to store preconditioner
  integer        :: l
!==============================================================================!

  !-----------------------------------------------------------!
  !   Fill up PETSc matrix with values from original matrix   !
  !-----------------------------------------------------------!
  do i = 1, Pet % m_lower
    do j = Sol % A % row(i), Sol % A % row(i+1)-1
      k = Sol % A % col(j)
      call C_Petsc_Mat_Set_Value(Pet % A,           &  ! matrix
                                 i-PETSC_ONE,       &  ! row
                                 k-PETSC_ONE,       &  ! column
                                 Sol % A % val(j))     ! value from orig. matrix
    end do
  end do

  ! The following two calls are needed after the calls to MatSetValue
  call C_Petsc_Assemble_Mat(Pet % A)

  !---------------------!
  !   Fill up vectors   !
  !---------------------!
  do i = 1, Pet % m_lower
    call C_Petsc_Vec_Set_Value(Pet % x, i - PETSC_ONE, x(i));
    call C_Petsc_Vec_Set_Value(Pet % b, i - PETSC_ONE, b(i));
  end do

  ! The following two calls are needed after the calls to VecSetValue
  call C_Petsc_Assemble_Vec(Pet % x)
  call C_Petsc_Assemble_Vec(Pet % b)

  !-----------------------------------!
  !   Set solver and preconditioner   !
  !-----------------------------------!
  sols  = KSPBICG;  l = len_trim(sols);   sols (l+1:l+1) = c_null_char
  precs = PCILU;    l = len_trim(precs);  precs(l+1:l+1) = c_null_char
  call C_Petsc_Set_Solver_And_Preconditioner(Pet % ksp,  &  ! solver
                                             Pet % pc,   &  ! preconditioner
                                             Pet % A,    &
                                             sols,       &
                                             precs)

  ! Set solver tolerances
  Pet % miter = miter
  call C_Petsc_Ksp_Set_Tolerances(Pet % ksp,     &
                                  tol,           &  ! PetscReal rtol
                                  tol,           &  ! PetscReal abstol
                                  Pet % miter)

  ! Solve
  call KSPSolve(Pet % ksp, Pet % b, Pet % x, Pet % err)

  ! Fetch the performed number of iterations
  call KSPGetIterationNumber(Pet % ksp, Pet % niter, Pet % err)
  niter = Pet % niter

  !-----------------------------------------------!
  !   Copy the solution back to T-Flows' vector   !
  !-----------------------------------------------!
  do i = 1, Pet % m_lower
    call C_Petsc_Vec_Get_Values(Pet % x,      &
                                PETSC_ONE,    &
                                i-PETSC_ONE,  &
                                x(i))
  end do

  end subroutine

  end module

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!                                                                              !
!   PETSc is not installed on the system -> take fake PETSc module             !
!                                                                              !
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
#else

!==============================================================================!
  module Petsc_Mod
!------------------------------------------------------------------------------!
!   This is a fake PETSc module, when code is compiled without it.             !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Solver_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !------------------!
  !   Solvers type   !
  !------------------!
  type Petsc_Type

    contains
      procedure :: Create_Petsc
      procedure :: Solve

  end type

  contains

!==============================================================================!
  subroutine Create_Petsc(Pet, Sol, Grid)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Petsc_Type)        :: Pet
  type(Solver_Type)        :: Sol
  type(Grid_Type),  target :: Grid
!==============================================================================!

  if(this_proc < 2) then
    print *, '# This version was compiled without PETSc,'
    print *, '# and yet they were specified in control file.'
    print *, '# This error is critical, exiting.'
  end if

  call Comm_Mod_End
  stop

  end subroutine

!==============================================================================!
  subroutine Solve(Pet, Sol, A, x, b, miter, niter, tol)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Petsc_Type)    :: Pet
  type(Solver_Type)    :: Sol
  type(Matrix_Type)    :: A
  real                 :: x(-Sol % pnt_grid % n_bnd_cells :  &
                             Sol % pnt_grid % n_cells)
  real                 :: b( Sol % pnt_grid % n_cells)
  integer, intent(in)  :: miter
  integer, intent(out) :: niter
  real,    intent(in)  :: tol
!==============================================================================!

  if(this_proc < 2) then
    print *, '# This version was compiled without PETSc,'
    print *, '# and yet they were specified in control file.'
    print *, '# This error is critical, exiting.'
  end if

  call Comm_Mod_End
  stop

  end subroutine

  end module

#endif
