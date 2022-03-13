!==============================================================================!
  subroutine Linear_Solvers(Rc, Flow, turb, Vof, Sol)
!------------------------------------------------------------------------------!
!   Reads which linear solvers to use, and calls appropriate branches.         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Read_Controls_Type)  :: Rc
  type(Field_Type),   target :: Flow
  type(Turb_Type),    target :: turb
  type(Vof_Type),     target :: Vof
  type(Solver_Type),  target :: Sol
!----------------------------------[Locals]------------------------------------!
  character(SL) :: name
!==============================================================================!

  ! Linear solvers you want to use; native or PETSc
  call Control_Mod_Linear_Solvers(name, .true.)
  Sol % solvers = Solver_Mod_Linear_Solvers_Code(name)

  ! Read options for native solvers first ...
  call Rc % Native_Solvers(Flow, turb, Vof, Sol)

  ! ... and follow with options for PETSc
  call Rc % Petsc_Solvers(Flow, turb, Vof, Sol)

  end subroutine
