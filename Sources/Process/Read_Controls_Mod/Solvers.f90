!==============================================================================!
  subroutine Solvers(Rc, Flow, Turb, Vof, Sol)
!------------------------------------------------------------------------------!
!   Reads which linear solvers to use, and calls appropriate branches.         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Read_Controls_Type), intent(in) :: Rc
  type(Field_Type),  target             :: Flow
  type(Turb_Type),   target             :: Turb
  type(Vof_Type),    target             :: Vof
  type(Solver_Type), target             :: Sol
!----------------------------------[Locals]------------------------------------!
  character(SL) :: name
!==============================================================================!

  ! Linear solvers you want to use; native or PETSc
  call Control % Linear_Solvers(name, .true.)
  Sol % solvers = Solver_Mod_Linear_Solvers_Code(name)

  if(Sol % solvers .eq. PETSC .and. .not. PETSC_ACTIVE) then
    call Message % Error(80,                                                 &
                  'This version was compiled without PETSc, and yet '    //  &
                  'they are specified in the control file.  This error ' //  &
                  'is critical, exiting. Either fix the control file  '  //  &
                  'by setting LINEAR_SOLVERS to native, or compile   '   //  &
                  'T-Flows with PETSc libraries.',                           &
                  file=__FILE__, line=__LINE__, one_proc=.true.)
  end if

  ! Read options for native solvers first ...
  call Rc % Native_Solvers(Flow, Turb, Vof)

  ! ... and follow with options for PETSc
  call Rc % Petsc_Solvers(Flow, Turb, Vof, Sol)

  end subroutine
