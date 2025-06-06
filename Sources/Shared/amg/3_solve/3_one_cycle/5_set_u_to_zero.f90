!==============================================================================!
  subroutine set_u_to_zero(amg, level, u)
!------------------------------------------------------------------------------!
!   Sets u-values of level "level" to zero
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(amg_type)  :: amg
  integer          :: level
  double precision :: u(:)
!-----------------------------------[locals]-----------------------------------!
  integer :: i
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  call amg % timer_start()

  do i = amg % imin(level), amg % imax(level)
    u(i) = 0.0d0
  end do

  call amg % timer_stop(15)

  end subroutine
