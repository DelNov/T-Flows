!==============================================================================!
  subroutine Calculate_Bulk_Velocities(Flow)
!------------------------------------------------------------------------------!
!   Calculate volume fluxes through whole domain.                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type), target :: Flow
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

  !$acc parallel loop
  do c = 1, grid_n_cells - grid_n_buff_cells
    bulk_u = bulk_u + u_n(c) * grid_vol(c)
    bulk_v = bulk_v + v_n(c) * grid_vol(c)
    bulk_w = bulk_w + w_n(c) * grid_vol(c)
    vol = vol + grid_vol(c)
  end do
  !$acc end parallel

  call Global % Sum_Real(bulk_u)
  call Global % Sum_Real(bulk_v)
  call Global % Sum_Real(bulk_w)
  call Global % Sum_Real(vol)

  ! Bulk velocities
  Flow % bulk % u = bulk_u / vol
  Flow % bulk % v = bulk_v / vol
  Flow % bulk % w = bulk_w / vol

  end subroutine
