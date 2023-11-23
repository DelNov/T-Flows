!==============================================================================!
  subroutine  Calculate_Global_Volumes(Grid)
!------------------------------------------------------------------------------!
!>  Calculates volume of the smallest and largest cell in a grid and total
!>  volume of the grid.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type) :: Grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c  ! counter
!==============================================================================!

  ! Initialize values
  Grid % min_vol =  HUGE
  Grid % max_vol = -HUGE
  Grid % tot_vol = 0.0

  ! Browse through cells avoiding buffers
  do c = 1, Grid % n_cells - Grid % Comm % n_buff_cells
    Grid % tot_vol = Grid % tot_vol + Grid % vol(c)
    Grid % min_vol = min(Grid % min_vol, Grid % vol(c))
    Grid % max_vol = max(Grid % max_vol, Grid % vol(c))
  end do

  ! Make them global (over all processors)
  call Global % Max_Real(Grid % max_vol)
  call Global % Min_Real(Grid % min_vol)
  call Global % Sum_Real(Grid % tot_vol)

  end subroutine
