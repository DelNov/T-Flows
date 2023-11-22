!==============================================================================!
  subroutine Remove_Singular(Sol, phi)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Solver_Type) :: Sol
  type(Var_Type)     :: phi
!==============================================================================!

  ! Remove singularity information for PETSc
  if(Sol % solvers == PETSC) then
    call C_Petsc_Mat_Remove_Null_Space(phi % Pet % A)
  end if

  end subroutine
