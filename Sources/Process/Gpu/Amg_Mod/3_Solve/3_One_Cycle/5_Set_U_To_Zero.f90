!==============================================================================!
  subroutine Set_U_To_Zero(Amg, level)
!------------------------------------------------------------------------------!
!   Sets u-values of level "level" to zero
!------------------------------------------------------------------------------!
  implicit none
!----------------------------------[Arguments]----------------------------------!
  class(Amg_Type),   target :: Amg
  integer                   :: level
!-----------------------------------[Locals]-----------------------------------!
  integer                   :: i, n
  real, contiguous, pointer :: u(:)
!------------------------------------[Save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  n =  Amg % lev(level) % n
  u => Amg % lev(level) % u

  do i = 1, n
    u(i) = 0.0
  end do

  end subroutine
