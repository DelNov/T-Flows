!==============================================================================!
  real function Cg_Alpha(Amg, level, s2,   &
                                     a, u, f, ia, ja,  &  ! defining system
                                     iw, m)
!------------------------------------------------------------------------------!
!   Called from Cg_Step, it seems to be updating alpha from the CG method
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(Amg_Type) :: Amg
  integer         :: level
  real            :: s2
  real            :: a(:), u(:), f(:)
  integer         :: ia(:), ja(:)
  integer         :: iw(:)
  integer         :: m
!-----------------------------------[locals]-----------------------------------!
  real    :: s1, sr
  integer :: i, iaux, ishift, j
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  ishift = Amg % imax(m) + 1 - Amg % imin(level)
  s1 = 0.0

  ! See comment in source "Coarsening.f90" at line 180
  iaux = ia(Amg % imax(level)+1)
  ia(Amg % imax(level)+1) = iw(Amg % iminw(level))

  do i = Amg % imin(level), Amg % imax(level)
    sr = 0.0
    do j = ia(i), ia(i+1) - 1
      sr = sr + a(j) * u(ja(j))
    end do
    s1 = s1 + sr * f(i+ishift)
  end do

  Cg_Alpha = -s1 / s2

  ia(Amg % imax(level)+1) = iaux

  end function
