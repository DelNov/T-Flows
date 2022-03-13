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

  logical, parameter :: PETSC_ACTIVE = .false.

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

  contains

!==============================================================================!
  subroutine Create_Petsc(Pet, Nat, Grid)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Petsc_Type)        :: Pet
  type(Native_Type)        :: Nat
  type(Grid_Type),  target :: Grid
!==============================================================================!

  Pet % pnt_grid => Grid

  if(this_proc < 2) then
    print *, '# NOTE! This version was compiled without PETSc, so the '  //  &
             'PETSc class was not created ...'
    print *, '# ... which is OK as long as you don''t specify to use ' //  &
             'PETSc solvers in control file.'
  end if

  end subroutine

!==============================================================================!
  subroutine Solve_Petsc(Pet,                      &
                         solver, prec, prec_opts,  &
                         A, x, b,                  &
                         miter, niter,             &
                         tol, fin_res,             &
                         norm)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Petsc_Type)    :: Pet
  character(*)         :: solver          ! solver
  character(*)         :: prec            ! preconditioner
  character(SL)        :: prec_opts(MSI)  ! preconditioner options
  type(Matrix_Type)    :: A
  real                 :: x(-Pet % pnt_grid % n_bnd_cells :  &
                             Pet % pnt_grid % n_cells)
  real                 :: b( Pet % pnt_grid % n_cells)
  integer, intent(in)  :: miter
  integer, intent(out) :: niter
  real,    intent(in)  :: tol      ! tolerance
  real,    intent(out) :: fin_res  ! final residual
  real,    optional    :: norm     ! normalization
!==============================================================================!

  if(this_proc < 2) then
    print *, '# This version was compiled without PETSc, '  //  &
             'and yet they were specified in the control file.'
    print *, '# This error is critical, exiting. Fix the control file.'
  end if

  call Comm_Mod_End
  stop

  end subroutine

  end module
