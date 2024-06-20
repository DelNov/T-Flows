!==============================================================================!
  subroutine Calculate_Bulk_Velocities(Flow, Grid)
!------------------------------------------------------------------------------!
!   Calculate volume fluxes through whole domain.                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type), target :: Flow
  type(Grid_Type)           :: Grid
!-----------------------------------[Locals]-----------------------------------!
  real    :: bulk_u, bulk_v, bulk_w, vol
  integer :: c
!==============================================================================!

  ! Set local variables which will not confuse OpenACC
  bulk_u = Flow % bulk % u
  bulk_v = Flow % bulk % v
  bulk_w = Flow % bulk % w

  ! Initialize bulk velocities to zero
  bulk_u = 0.0
  bulk_v = 0.0
  bulk_w = 0.0

  !----------------------------------------------------------------------!
  !   Summ up volume fluxes [m^3/s] over all faces at monitoring plane   !
  !----------------------------------------------------------------------!

  vol = 0.0

  !$acc parallel loop reduction(+:bulk_u, bulk_v, bulk_w, vol)  &
  !$acc present(grid_region_f_cell, grid_region_l_cell,         &
  !$acc         grid_vol,                                       &
  !$acc         flow_u_n, flow_v_n, flow_w_n)
  do c = Cells_In_Domain_Gpu()  ! all present
    bulk_u = bulk_u + flow_u_n(c) * grid_vol(c)
    bulk_v = bulk_v + flow_v_n(c) * grid_vol(c)
    bulk_w = bulk_w + flow_w_n(c) * grid_vol(c)
    vol = vol + grid_vol(c)
  end do
  !$acc end parallel

  call Global % Sum_Reals(bulk_u, bulk_v, bulk_w, vol)

  ! Bulk velocities
  Flow % bulk % u = bulk_u / vol
  Flow % bulk % v = bulk_v / vol
  Flow % bulk % w = bulk_w / vol

  end subroutine
