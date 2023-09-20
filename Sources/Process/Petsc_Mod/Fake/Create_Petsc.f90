!==============================================================================!
  subroutine Create_Petsc(Pet, A, var_name, options_pets)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Petsc_Type) :: Pet
  type(Matrix_Type) :: A
  character(VL)     :: var_name
  character(SL)     :: options_pets(MAX_STRING_ITEMS)
!-----------------------------------[Locals]-----------------------------------!
  logical, save :: called = .false.
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

