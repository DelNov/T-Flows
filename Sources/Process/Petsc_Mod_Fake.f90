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
