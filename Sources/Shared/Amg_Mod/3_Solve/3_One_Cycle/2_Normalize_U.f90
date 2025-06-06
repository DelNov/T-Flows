!==============================================================================!
  subroutine Normalize_U(Amg, level, u)
!------------------------------------------------------------------------------!
!   Normalizes u on level "level" if rowsum=0 (last component =0)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(Amg_Type)  :: Amg
  integer          :: level
  double precision :: u(:)
!-----------------------------------[locals]-----------------------------------!
  double precision :: fac
  integer          :: i
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  if(Amg % irow0 .eq. AMG_NON_SINGULAR_MATRIX) return

  fac = u(Amg % imax(level))

  do i = Amg % imin(level), Amg % imax(level)
    u(i) = u(i) - fac
  end do

  end subroutine
