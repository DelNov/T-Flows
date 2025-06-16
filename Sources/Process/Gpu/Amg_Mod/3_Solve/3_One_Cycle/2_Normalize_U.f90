!==============================================================================!
  subroutine Normalize_U(Amg, level)
!------------------------------------------------------------------------------!
!   Normalizes u on level "level" if rowsum=0 (last component =0)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(Amg_Type), target :: Amg
  integer                 :: level
!-----------------------------------[locals]-----------------------------------!
  real                      :: fac
  integer                   :: i, n
  real, contiguous, pointer :: u(:)
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  if(Amg % irow0 .eq. AMG_NON_SINGULAR_MATRIX) return

  n  =  Amg % lev(level) % n
  u  => Amg % lev(level) % u

  fac = u(n)

  do i = 1, n
    u(i) = u(i) - fac
  end do

  end subroutine
