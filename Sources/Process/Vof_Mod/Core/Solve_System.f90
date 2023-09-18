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
!==============================================================================!

  ! Take aliases
  fun  => Vof % fun
  Flow => Vof % pnt_flow
  A    => Sol % Nat % A

  call Profiler % Start(String % First_Upper(fun % solver)//' (solver for VOF)')

  ! Call linear solver to solve the equations
  call Sol % Run(A, fun, b)

  call Profiler % Stop(String % First_Upper(fun % solver)//' (solver for VOF)')

  if(.not. Flow % heat_transfer) then
    call Info % Iter_Fill_At(1, 6, fun % name, fun % res, fun % niter)
  else
    call Info % Iter_Fill_At(2, 1, fun % name, fun % res, fun % niter)
  end if

  end subroutine
