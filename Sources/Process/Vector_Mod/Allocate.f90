!==============================================================================!
  subroutine Vector_Mod_Allocate(grid, vector)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)   :: grid
  type(Vector_Type) :: vector
!==============================================================================!

  ! Allocate memory for vector
  allocate (vector % val(grid % n_cells));  vector % val = 0.0

  end subroutine
