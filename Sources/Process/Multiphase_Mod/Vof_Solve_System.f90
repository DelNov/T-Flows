!==============================================================================!
  subroutine Multiphase_Mod_Vof_Solve_System(mult, sol, b)
!------------------------------------------------------------------------------!
!   Solves linear system for VOF                                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Solver_Type),     target :: sol
  type(Multiphase_Type), target :: mult
  real,                  target :: b(:)
!-----------------------------------[Locals]-----------------------------------!
  type(Var_Type),   pointer :: vof
  character(len=80)         :: solver
  integer                   :: exec_iter
!==============================================================================!

  ! Take aliases
  vof => mult % vof

  ! Get solver
  call Control_Mod_Solver_For_Multiphase(solver)

  call Cpu_Timer_Mod_Start('Linear_Solver_For_Multiphase')
  call Cg(sol,            &
          vof % n,        &
          b,              &
          vof % precond,  &
          vof % niter,    &      ! max number of iterations
          exec_iter,      &      ! executed number of iterations
          vof % tol,      &
          vof % res)
  call Cpu_Timer_Mod_Stop('Linear_Solver_For_Multiphase')

  call Info_Mod_Iter_Fill_At(1, 6, vof % name, exec_iter, vof % res)

  end subroutine
