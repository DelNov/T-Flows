!==============================================================================!
  subroutine Create_Solver(Sol, Grid)
!------------------------------------------------------------------------------!
!>  The Create_Solver subroutine within the Solver_Mod module is designed for
!>  initializing the solver environment within T-Flows. This subroutine is only
!>  focused on setting up T-Flows' native solver, as the PETSc solvers are
!>  encapsulated and managed separately within the Var_Mod module, where they
!>  are initialized in conjunction with variable creation. This routine's
!>  exclusive focus on native solvers (in contrast to the End subroutine's
!>  focus on PETSc solvers) reflects the non-symmeric nature of the Solver_Mod
!>  module in handling both native and PETSc solvers within T-Flows.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Solver_Type)       :: Sol   !! parent class
  type(Grid_Type),  target :: Grid  !! grid on which the solver will be used
!==============================================================================!

  call Sol % Nat % Create_Native(Grid)

  end subroutine

