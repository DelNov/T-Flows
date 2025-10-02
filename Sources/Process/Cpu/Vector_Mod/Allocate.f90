!==============================================================================!
  subroutine Vector_Mod_Allocate(vector, Grid)
!------------------------------------------------------------------------------!
!>  This subroutine is responsible for allocating memory and initializing the
!>  vector data structure within T-Flows. The size of the vector is set to be
!>  equal to the number of computational cells in the grid and the link to the
!>  computational grid is established.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Vector_Type)        :: vector  !! vector being allocated
  type(Grid_Type),  target :: Grid    !! computational grid on which
                                      !! the vector is being defined
!==============================================================================!

  ! Store pointer to the grid
  vector % pnt_grid => Grid

  ! Allocate memory for vector
  allocate (vector % val(Grid % n_cells));  vector % val = 0.0

  end subroutine
