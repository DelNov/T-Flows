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

  ! First condition is to see if you compiled with PETSc at all,
  ! the second is to check if PETSc was engaged by the user

# if T_FLOWS_PETSC == 1
    if(PETSC_ACTIVE) then
      if(petsc_is_reporting) call C_Petsc_Log_View()
      call C_Petsc_Finalize()
    end if
# endif

  end subroutine
