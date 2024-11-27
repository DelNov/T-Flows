!==============================================================================!
  subroutine Destroy_Petsc(Pet, var_name)
!------------------------------------------------------------------------------!
!>  This subroutine serves as placeholders or dummy variant when PETSc is
!>  not available in the system.  The functionality here is intentionally
!>  minimal, primarily serving to maintain code structure and compatibility
!>  in environments where PETSc is not used.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Petsc_Type) :: Pet       !! parent object of the Petsc_Type
  character(VL)     :: var_name  !! variable name
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Pet)
  Unused(var_name)
!==============================================================================!

  end subroutine

