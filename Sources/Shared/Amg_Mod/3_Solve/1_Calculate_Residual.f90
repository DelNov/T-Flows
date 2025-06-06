!==============================================================================!
  subroutine Calculate_Residual(Amg, level, resl,  &
                                a, u, f, ia, ja,   &  ! defining system
                                iw)
!------------------------------------------------------------------------------!
!   Computes l2-norm of residual on level "level"
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(Amg_Type)  :: Amg
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

  ! See comment in source "Coarsening.f90" at line 180
  iaux = ia(Amg % imax(level)+1)
  ia(Amg % imax(level)+1) = iw(Amg % iminw(level))

  do i = Amg % imin(level), Amg % imax(level)
    s = f(i)
    do j = ia(i), ia(i+1) - 1
      s = s - a(j) * u(ja(j))
    end do
    resl = resl + s*s
  end do

  ia(Amg % imax(level)+1) = iaux
  resl = sqrt(resl)

  end subroutine
