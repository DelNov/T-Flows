!==============================================================================!
  subroutine Vector_Mod_Allocate_Level(vector, grid, lev)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Vector_Type)       :: vector
  type(Grid_Type), target :: grid
  integer                 :: lev
!==============================================================================!

  ! Store pointer to the grid
  vector % pnt_grid => grid

  ! Allocate memory for vector
  allocate (vector % val(grid % level(lev) % n_cells));  vector % val(:) = 0.0

  end subroutine
