!==============================================================================!
  subroutine Create_Solver(Sol, grid)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Solver_Type)       :: Sol
  type(Grid_Type),  target :: grid
!==============================================================================!

  Sol % pnt_grid => grid

  if(this_proc < 2) print *, '# Determining matrix topology.'

  call Sol % A % Create_Matrix(grid)
  call Sol % M % Create_Matrix(grid)
  call Sol % D % Create_Matrix(grid)
  call Vector_Mod_Allocate(Sol % b, grid)

  if(this_proc < 2) print *, '# Finished !'

  end subroutine
