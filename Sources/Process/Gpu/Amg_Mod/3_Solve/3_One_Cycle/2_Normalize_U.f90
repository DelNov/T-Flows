!==============================================================================!
  subroutine Normalize_U(Amg, level)
!------------------------------------------------------------------------------!
!   Normalizes u on level "level" if rowsum=0 (last component =0)
!------------------------------------------------------------------------------!
  implicit none
!----------------------------------[Arguments]----------------------------------!
  class(Amg_Type), target :: Amg
  integer                 :: level
!-----------------------------------[Locals]-----------------------------------!
  real                      :: fac
  integer                   :: i, n
  real, contiguous, pointer :: u(:)
!------------------------------------[Save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  if(Amg % matrix % singular .eq. AMG_NON_SINGULAR_MATRIX) return

  n  =  Amg % lev(level) % n
  u  => Amg % lev(level) % u

  fac = u(n)

  do i = 1, n
    u(i) = u(i) - fac
  end do

  end subroutine
