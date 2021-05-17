!==============================================================================!
  subroutine Create_Solver(Sol, Grid)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Solver_Type)       :: Sol
  type(Grid_Type),  target :: Grid
!==============================================================================!

  Sol % pnt_grid => Grid

  if(this_proc < 2) print *, '# Determining matrix topology.'

  call Sol % A % Create_Matrix(Grid)
  call Sol % M % Create_Matrix(Grid)
  call Sol % D % Create_Matrix(Grid)
  call Vector_Mod_Allocate(Sol % b, Grid)

  if(this_proc < 2) print *, '# Finished !'

  end subroutine
