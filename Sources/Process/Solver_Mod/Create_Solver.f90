!==============================================================================!
  subroutine Create_Solver(Sol, Grid)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Solver_Type)       :: Sol
  type(Grid_Type),  target :: Grid
!==============================================================================!

  call Sol % Nat % Create_Native(Grid)
  call Sol % Pet % Create_Petsc(Sol % Nat, Grid)

  end subroutine

