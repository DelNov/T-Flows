!==============================================================================!
  subroutine normalize_u(amg, level, u)
!------------------------------------------------------------------------------!
!   Normalizes u on level "level" if rowsum=0 (last component =0)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(amg_type)  :: amg
  integer          :: level
  double precision :: u(:)
!-----------------------------------[locals]-----------------------------------!
  double precision :: fac
  integer          :: i
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  if(amg % irow0 .eq. AMG_NON_SINGULAR_MATRIX) return

  fac = u(amg % imax(level))

  do i = amg % imin(level), amg % imax(level)
    u(i) = u(i) - fac
  end do

  end subroutine
