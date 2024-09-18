!==============================================================================!
  real function Volume_Average(Flow, Grid, val)
!------------------------------------------------------------------------------!
!   Calculates the volume average of the values in array val                   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type), target :: Flow  !! parent flow field class
  type(Grid_Type),   target :: Grid  !! grid on which the flow is defined
  real                      :: val(-Grid % n_bnd_cells:Grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  integer :: c
  real    :: sum_val, tot_vol
!==============================================================================!

  ! Initialize sums
  sum_val = 0.0
  tot_vol = 0.0

  !-------------------------------------------------!
  !   Integrate them over cells, avoiding buffers   !
  !-------------------------------------------------!

  ! For this to work properly on GPUs, one should make sure that variable
  ! "val" is coppied to GPU.  So far it was only density which was sent
  ! here, so things should work.  But it is probably good to stay alert.

  !$acc parallel loop independent reduction(+: sum_val,tot_vol)  &
  !$acc present(  &
  !$acc   grid_region_f_cell,  &
  !$acc   grid_region_l_cell,  &
  !$acc   val,  &
  !$acc   grid_vol   &
  !$acc )
  do c = grid_region_f_cell(grid_n_regions), grid_region_l_cell(grid_n_regions)  ! all present
    sum_val = sum_val + val(c) * grid_vol(c)
    tot_vol = tot_vol + grid_vol(c)
  end do
  !$acc end parallel

  ! Perform sums over all processors
  call Global % Sum_Real(sum_val)
  call Global % Sum_Real(tot_vol)

  ! Work out the final resul
  Volume_Average = sum_val / tot_vol

  end function
