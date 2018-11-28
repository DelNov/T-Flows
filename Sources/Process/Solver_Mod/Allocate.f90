!==============================================================================!
  subroutine Solver_Mod_Allocate(sol, grid)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Solver_Type)        :: sol
  type(Grid_Type),  target :: grid
!==============================================================================!

  sol % pnt_grid => grid

  call Matrix_Mod_Create(grid, sol % d)
  call Matrix_Mod_Create(grid, sol % a)
  allocate (sol % b(grid % n_cells));    sol % b(:) = 0.0

  end subroutine

