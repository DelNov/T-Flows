!==============================================================================!
  subroutine Remove_Singular(Sol, A)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Solver_Type) :: Sol
  type(Matrix_Type)  :: A
!==============================================================================!

  ! Remove singularity information for PETSc
  if(Sol % solvers == PETSC) then
    call C_Petsc_Mat_Remove_Null_Space(Sol % Pet % A)
  end if

  end subroutine
