!==============================================================================!
  subroutine Create_Petsc(Pet, A, var_name, petsc_rank)
!------------------------------------------------------------------------------!
!>  This subroutine serves as placeholders or dummy variant when PETSc is
!>  not available in the system.  The functionality here is intentionally
!>  minimal, primarily serving to maintain code structure and compatibility
!>  in environments where PETSc is not used.  The only functionality in this
!>  subroutine is to issue a warning if a user tries to use PETSc solvers (by
!>  specifying so in the control file), but they are not available in the code.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Petsc_Type)    :: Pet         !! parent object of the Petsc_Type
  type(Matrix_Type)    :: A           !! matrix in T-Flows native format
  character(VL)        :: var_name    !! variable name
  integer, intent(out) :: petsc_rank  !! rank of the Petsc object instance
!-----------------------------------[Locals]-----------------------------------!
  logical, save :: called = .false.
!------------------------------[Unused arguments]------------------------------!
  Unused(var_name)
  Unused(petsc_rank)
!==============================================================================!

  Pet % pnt_grid => A % pnt_grid

  if(.not. called) then
    if(First_Proc()) then
      print '(a)', ' # NOTE! This version was compiled without PETSc,'  //  &
                   ' so the PETSc class was not created ...'
      print '(a)', ' # ... which is OK as long as you don''t specify'   //  &
                   ' to use PETSc solvers in control file.'
    end if
    called = .true.
  end if

  end subroutine

