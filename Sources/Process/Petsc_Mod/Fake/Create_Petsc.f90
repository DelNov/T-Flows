!==============================================================================!
  subroutine Create_Petsc(Pet, Nat, Grid)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Petsc_Type)        :: Pet
  type(Native_Type)        :: Nat
  type(Grid_Type),  target :: Grid
!==============================================================================!

  Pet % pnt_grid => Grid

  if(this_proc < 2) then
    print '(a)', ' # NOTE! This version was compiled without PETSc,'  //  &
                 ' so the PETSc class was not created ...'
    print '(a)', ' # ... which is OK as long as you don''t specify'   //  &
                 ' to use PETSc solvers in control file.'
  end if

  end subroutine

