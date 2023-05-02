!==============================================================================!
  subroutine End(Sol)
!------------------------------------------------------------------------------!
!   End PETSc.  If this is not called, an error is thrown on some systems.     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Solver_Type) :: Sol
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Sol)
!==============================================================================!

  ! Call linear solver to solve the equations
  if(PETSC_ACTIVE) then
    call C_Petsc_Finalize()
  end if

  end subroutine
