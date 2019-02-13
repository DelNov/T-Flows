!==============================================================================!
  subroutine Vector_Mod_Allocate(vector, grid)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Vector_Type)        :: vector
  type(Grid_Type),  target :: grid
!==============================================================================!

  ! Store pointer to the grid
  vector % pnt_grid => grid

  ! Allocate memory for vector
  allocate (vector % val(grid % n_cells));  vector % val = 0.0

  end subroutine
