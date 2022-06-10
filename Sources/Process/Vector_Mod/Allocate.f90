!==============================================================================!
  subroutine Vector_Mod_Allocate(vector, Grid)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Vector_Type)        :: vector
  type(Grid_Type),  target :: Grid
!==============================================================================!

  ! Store pointer to the grid
  vector % pnt_grid => Grid

  ! Allocate memory for vector
  allocate (vector % val(Grid % n_cells));  vector % val = 0.0

  end subroutine
