!==============================================================================!
  subroutine Multiphase_Mod_Vof_Solve_System(mult, Sol, b)
!------------------------------------------------------------------------------!
!   Solves linear system for VOF                                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Solver_Type),     target :: Sol
  type(Multiphase_Type), target :: mult
  real, contiguous,      target :: b(:)
!-----------------------------------[Locals]-----------------------------------!
  type(Var_Type),    pointer :: vof
  type(Field_Type),  pointer :: flow
  type(Matrix_Type), pointer :: A
  character(SL)              :: solver
!==============================================================================!

  ! Take aliases
  vof  => mult % vof
  flow => mult % pnt_flow
  A    => Sol % A

  ! Get solver
  call Control_Mod_Solver_For_Multiphase(solver)

  call Cpu_Timer % Start('Linear_Solver_For_Multiphase')
  call Sol % Bicg(A,              &
                  vof % n,        &
                  b,              &
                  vof % precond,  &
                  vof % mniter,   &      ! max number of iterations
                  vof % eniter,   &      ! executed number of iterations
                  vof % tol,      &
                  vof % res)
  call Cpu_Timer % Stop('Linear_Solver_For_Multiphase')

  if(.not. flow % heat_transfer) then
    call Info_Mod_Iter_Fill_At(1, 6, vof % name, vof % eniter, vof % res)
  else
    call Info_Mod_Iter_Fill_At(2, 1, vof % name, vof % eniter, vof % res)
  end if

  end subroutine
