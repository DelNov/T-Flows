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

  !$acc parallel loop
  do c = Cells_In_Domain()
    bulk_u = bulk_u + u_n(c) * Grid % vol(c)
    bulk_v = bulk_v + v_n(c) * Grid % vol(c)
    bulk_w = bulk_w + w_n(c) * Grid % vol(c)
    vol = vol + Grid % vol(c)
  end do
  !$acc end parallel

  call Global % Sum_Reals(bulk_u, bulk_v, bulk_w, vol)

  ! Bulk velocities
  Flow % bulk % u = bulk_u / vol
  Flow % bulk % v = bulk_v / vol
  Flow % bulk % w = bulk_w / vol

  end subroutine