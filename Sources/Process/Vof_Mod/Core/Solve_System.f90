!==============================================================================!
  subroutine Solve_System(Vof, Sol, b)
!------------------------------------------------------------------------------!
!   Solves linear system for VOF                                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type),   target :: Vof
  type(Solver_Type), target :: Sol
  real, contiguous,  target :: b(:)
!-----------------------------------[Locals]-----------------------------------!
  type(Var_Type),    pointer :: fun
  type(Field_Type),  pointer :: Flow
  type(Matrix_Type), pointer :: A
  character(SL)              :: solver
!==============================================================================!

  ! Take aliases
  fun  => Vof % fun
  Flow => Vof % pnt_flow
  A    => Sol % Nat % A

  ! Get solver
  call Control_Mod_Solver_For_Vof(solver)

  call Profiler % Start('Linear_Solver_For_Vof')

  ! Call linear solver to solve the equations
  call Sol % Run(fun % solver,     &
                 fun % prec,       &
                 fun % prec_opts,  &
                 A,                &
                 fun % n,          &
                 b,                &
                 fun % mniter,     &
                 fun % eniter,     &
                 fun % tol,        &
                 fun % res)

  call Profiler % Stop('Linear_Solver_For_Vof')

  if(.not. Flow % heat_transfer) then
    call Info_Mod_Iter_Fill_At(1, 6, fun % name, fun % eniter, fun % res)
  else
    call Info_Mod_Iter_Fill_At(2, 1, fun % name, fun % eniter, fun % res)
  end if

  end subroutine
