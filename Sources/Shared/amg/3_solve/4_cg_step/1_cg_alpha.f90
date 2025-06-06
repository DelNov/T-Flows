!==============================================================================!
  double precision function cg_alpha(amg, level, s2,   &
                                     a, u, f, ia, ja,  &  ! defining system
                                     iw, m)
!------------------------------------------------------------------------------!
!   Called from cg_step, it seems to be updating alpha from the CG method
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(amg_type)  :: amg
  integer          :: level
  double precision :: s2
  double precision :: a(:), u(:), f(:)
  integer          :: ia(:), ja(:)
  integer          :: iw(:)
  integer          :: m
!-----------------------------------[locals]-----------------------------------!
  double precision :: s1, sr
  integer          :: i, iaux, ishift, j
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  ishift = amg % imax(m) + 1 - amg % imin(level)
  s1 = 0.0d0

  ! See comment in source "coarsening.f90" at line 180
  iaux = ia(amg % imax(level)+1)
  ia(amg % imax(level)+1) = iw(amg % iminw(level))

  do i = amg % imin(level), amg % imax(level)
    sr = 0.0d0
    do j = ia(i), ia(i+1) - 1
      sr = sr + a(j) * u(ja(j))
    end do
    s1 = s1 + sr * f(i+ishift)
  end do

  cg_alpha = -s1 / s2

  ia(amg % imax(level)+1) = iaux

  end function
