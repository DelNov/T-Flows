!==============================================================================!
  subroutine Backup_U(Amg, level, icgr)
!------------------------------------------------------------------------------!
!   Makes a back-up of the current approx. on level "level" if icgr.ne.0.
!   (This seems to try to place backup beyond the last (coarsest) level.)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(Amg_Type), target :: Amg
  integer                 :: level, icgr
!-----------------------------------[locals]-----------------------------------!
  integer                      :: i, n
  real,    contiguous, pointer :: u(:), u_b(:)
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  if(icgr .eq. 0) return

  n   =  Amg % lev(level) % n
  u   => Amg % lev(level) % u
  u_b => Amg % lev(level) % u_b

  do i = 1, n
    u_b(i) = u(i)
  end do

  end subroutine
