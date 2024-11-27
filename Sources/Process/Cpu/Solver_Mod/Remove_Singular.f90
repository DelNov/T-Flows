!==============================================================================!
  subroutine Remove_Singular(Sol, phi)
!------------------------------------------------------------------------------!
!>  The Remove_Singular subroutine's main focus is on managing singularities
!>  in matrices when using PETSc solvers.  Its main purpose is to remove
!>  singularity from the matrix associated with a specific variable (phi),
!>  but only when the PETSc solver is in use.  For native solvers removing
!>  singularity wouldn't make much sense, since there is only instance of
!>  native solvers (and its matrices) and matrix is re-formed after the
!>  solution anyhow.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Solver_Type) :: Sol  !! parent class
  type(Var_Type)     :: phi  !! unknown variable
!==============================================================================!

  ! Remove singularity information for PETSc
  if(Sol % solvers == PETSC) then
    call C_Petsc_Mat_Remove_Null_Space(phi % Pet % A)
  end if

  end subroutine
