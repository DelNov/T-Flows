!==============================================================================!
  subroutine Allocate_Nodes(Grid, nn)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type) :: Grid
  integer          :: nn
!-----------------------------------[Locals]-----------------------------------!
  integer :: n
!==============================================================================!

  ! Store the number of nodes for the grid
  Grid % n_nodes = nn

  ! Allocate memory for node coordinates
  allocate(Grid % xn(1:nn));  Grid % xn(:) = 0.0
  allocate(Grid % yn(1:nn));  Grid % yn(:) = 0.0
  allocate(Grid % zn(1:nn));  Grid % zn(:) = 0.0

  allocate(Grid % Comm % node_glo(1:nn));  Grid % Comm % node_glo(:) = 0
  do n = 1, nn
    Grid % Comm % node_glo(n) = n
  end do

  ! Allocate new and old numbers (this is so often used, maybe is better here)
  allocate(Grid % new_n(n));  Grid % new_n(:) = 0
  allocate(Grid % old_n(n));  Grid % old_n(:) = 0

  end subroutine
