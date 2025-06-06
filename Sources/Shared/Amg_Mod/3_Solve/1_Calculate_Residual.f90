!==============================================================================!
  subroutine calculate_residual(amg, level, resl,  &
                                a, u, f, ia, ja,   &  ! defining system
                                iw)
!------------------------------------------------------------------------------!
!   Computes l2-norm of residual on level "level"
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(amg_type)  :: amg
  integer          :: level
  double precision :: resl
  double precision :: a(:), u(:), f(:)
  integer          :: ia(:), ja(:)
  integer          :: iw(:)
!-----------------------------------[locals]-----------------------------------!
  double precision :: s
  integer          :: i, iaux, j
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  resl = 0.0d0

  ! See comment in source "coarsening.f90" at line 180
  iaux = ia(amg % imax(level)+1)
  ia(amg % imax(level)+1) = iw(amg % iminw(level))

  do i = amg % imin(level), amg % imax(level)
    s = f(i)
    do j = ia(i), ia(i+1) - 1
      s = s - a(j) * u(ja(j))
    end do
    resl = resl + s*s
  end do

  ia(amg % imax(level)+1) = iaux
  resl = sqrt(resl)

  end subroutine
