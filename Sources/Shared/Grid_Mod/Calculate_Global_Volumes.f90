!==============================================================================!
  subroutine  Calculate_Global_Volumes(Grid)
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
  do c = 1, Grid % n_cells - Grid % comm % n_buff_cells
    Grid % tot_vol = Grid % tot_vol + Grid % vol(c)
    Grid % min_vol = min(Grid % min_vol, Grid % vol(c))
    Grid % max_vol = max(Grid % max_vol, Grid % vol(c))
  end do

  ! Make them global (over all processors)
  call Comm_Mod_Global_Max_Real(Grid % max_vol)
  call Comm_Mod_Global_Min_Real(Grid % min_vol)
  call Comm_Mod_Global_Sum_Real(Grid % tot_vol)

  end subroutine
