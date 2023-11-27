!==============================================================================!
  subroutine End(Sol)
!------------------------------------------------------------------------------!
!>  This subroutine is tailored for correctly terminating the PETSc environment.
!>  Calling End at the appropriate time in the simulation lifecycle is crucial
!>  for ensuring a clean and error-free shutdown of PETSc. This routine's
!>  exclusive focus on PETSc (in contrast to the Create_Solver subroutine's
!>  focus on native solvers) reflects the non-symmeric nature of the Solver_Mod
!>  module in handling both native and PETSc solvers within T-Flows.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Solver_Type) :: Sol  !! parent class
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Sol)
!==============================================================================!

  ! First condition is to see if you compiled with PETSc at all,
  ! the second is to check if PETSc was engaged by the user

# if T_FLOWS_PETSC == 1
    if(PETSC_ACTIVE) then
      if(petsc_is_reporting) call C_Petsc_Log_View()
      call C_Petsc_Finalize()
    end if
# endif

  end subroutine
