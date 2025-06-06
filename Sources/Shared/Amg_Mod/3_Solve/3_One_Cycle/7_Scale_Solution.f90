!==============================================================================!
  subroutine Scale_Solution(Amg, level, ivstar,  &
                            a, u, f, ia, ja,     &  ! defining system
                            iw)
!------------------------------------------------------------------------------!
!   Scales actual approximate solution on level "level" (v*-cycle); scaling
!   is done such that energy norm becomes minimal
!
!   Note: this scaling makes sense only for symmetric problems
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(Amg_Type)  :: Amg
  integer          :: level, ivstar
  double precision :: a(:), u(:), f(:)
  integer          :: ia(:), ja(:)
  integer          :: iw(:)
!-----------------------------------[locals]-----------------------------------!
  double precision :: fac, s1, s2, sa
  integer          :: i, iaux, j
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  if(ivstar .ne. 1) return

  call Amg % timer_start()

  !-----------------------------------------!
  !   Computation of scaling factor "fac"   !
  !-----------------------------------------!

  ! See comment in source "Coarsening.f90" at line 180
  iaux = ia(Amg % imax(level)+1)
  ia(Amg % imax(level)+1) = iw(Amg % iminw(level))

  s1 = 0.0d0
  s2 = 0.0d0
  do i = Amg % imin(level), Amg % imax(level)
    sa = 0.0d0
    do j = ia(i), ia(i+1) - 1
      sa = sa+a(j)*u(ja(j))
    end do
    s1 = s1+u(i)*f(i)
    s2 = s2+u(i)*sa
  end do

  fac = 1.0d0
  if (s2.ne.0.0d0) fac = s1/s2

  !-------------!
  !   Scaling   !
  !-------------!
  do i = Amg % imin(level), Amg % imax(level)
    u(i) = u(i)*fac
  end do
  ia(Amg % imax(level)+1) = iaux

  call Amg % timer_stop(14)

  end subroutine
