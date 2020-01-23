!==============================================================================!
  subroutine  Grid_Mod_Calculate_Global_Volumes(grid)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c  ! counter
!==============================================================================!

  ! Initialize values
  grid % min_vol =  HUGE
  grid % max_vol = -HUGE
  grid % tot_vol = 0.0

  ! Browse through cells avoiding buffers
  do c = 1, grid % n_cells - grid % comm % n_buff_cells
    grid % tot_vol = grid % tot_vol + grid % vol(c)
    grid % min_vol = min(grid % min_vol, grid % vol(c))
    grid % max_vol = max(grid % max_vol, grid % vol(c))
  end do

  ! Make them global (over all processors)
  call Comm_Mod_Global_Max_Real(grid % max_vol)
  call Comm_Mod_Global_Min_Real(grid % min_vol)
  call Comm_Mod_Global_Sum_Real(grid % tot_vol)

  end subroutine
