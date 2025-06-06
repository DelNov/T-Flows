!==============================================================================!
  subroutine Set_U_To_Zero(Amg, level, u)
!------------------------------------------------------------------------------!
!   Sets u-values of level "level" to zero
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(Amg_Type)  :: Amg
  integer          :: level
  double precision :: u(:)
!-----------------------------------[locals]-----------------------------------!
  integer :: i
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  call Amg % timer_start()

  do i = Amg % imin(level), Amg % imax(level)
    u(i) = 0.0d0
  end do

  call Amg % timer_stop(15)

  end subroutine
