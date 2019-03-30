!==============================================================================!
  subroutine Grid_Mod_Allocate_Nodes(grid, n)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: n
!==============================================================================!

  ! Store the number of nodes for the grid
  grid % n_nodes = n

  ! Allocate memory for node coordinates
  allocate(grid % xn(1:n));  grid % xn = 0.0
  allocate(grid % yn(1:n));  grid % yn = 0.0
  allocate(grid % zn(1:n));  grid % zn = 0.0

  end subroutine
