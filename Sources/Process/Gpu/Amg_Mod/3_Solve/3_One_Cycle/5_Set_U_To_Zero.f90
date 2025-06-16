!==============================================================================!
  subroutine Set_U_To_Zero(Amg, level, u)
!------------------------------------------------------------------------------!
!   Sets u-values of level "level" to zero
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(Amg_Type),   target :: Amg
  integer                   :: level
  real                      :: u(:)
  real, contiguous, pointer :: lev_u(:)
!-----------------------------------[locals]-----------------------------------!
  integer :: i, n
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  call Amg % timer_start()

!old:  do i = Amg % imin(level), Amg % imax(level)
!old:    u(i) = 0.0
!old:  end do

  n      =  Amg % lev(level) % n
  lev_u  => Amg % lev(level) % u

  ! Copy vectors u and f to level's storage
  ! (This is needed during the development stage)
  call Amg % Update_U_And_F_At_Level(level, vec_u=u)

  do i = 1, n
    lev_u(i) = 0.0
  end do

  ! Copy vectors u and f back to global storage
  ! (This is needed during the development stage)
  call Amg % Update_U_And_F_Globally(level, vec_u=u)

  call Amg % timer_stop(15)

  end subroutine
