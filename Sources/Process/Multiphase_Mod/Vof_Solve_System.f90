!==============================================================================!
  subroutine Multiphase_Mod_Vof_Solve_System(mult, sol, b)
!------------------------------------------------------------------------------!
!   Solves linear system for VOF                                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Solver_Type),     target :: sol
  type(Multiphase_Type), target :: mult
  real, contiguous,      target :: b(:)
!-----------------------------------[Locals]-----------------------------------!
  type(Var_Type),    pointer :: vof
  type(Matrix_Type), pointer :: a
  character(SL)              :: solver
!==============================================================================!

  ! Take aliases
  vof => mult % vof
  a   => sol % a

  ! Get solver
  call Control_Mod_Solver_For_Multiphase(solver)

  call Cpu_Timer_Mod_Start('Linear_Solver_For_Multiphase')
  call Solver_Mod_Bicg(sol,            &
                       vof % n,        &
                       b,              &
                       vof % precond,  &
                       vof % mniter,   &      ! max number of iterations
                       vof % eniter,   &      ! executed number of iterations
                       vof % tol,      &
                       vof % res)
  call Cpu_Timer_Mod_Stop('Linear_Solver_For_Multiphase')

  if(.not. mult % mass_transfer) then
    call Info_Mod_Iter_Fill_At(1, 6, vof % name, vof % eniter, vof % res)
  end if

  end subroutine
