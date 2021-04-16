!==============================================================================!
  subroutine Solver_Mod_Create(sol, grid)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Solver_Type)        :: sol
  type(Grid_Type),  target :: grid
!==============================================================================!

  sol % pnt_grid => grid

  if(this_proc < 2) print *, '# Determining matrix topology.'

  call Matrix_Mod_Create  (sol % a, grid)
  call Matrix_Mod_Create  (sol % m, grid)
  call Matrix_Mod_Create  (sol % d, grid)
  call Vector_Mod_Allocate(sol % b, grid)

  if(this_proc < 2) print *, '# Finished !'

  end subroutine
