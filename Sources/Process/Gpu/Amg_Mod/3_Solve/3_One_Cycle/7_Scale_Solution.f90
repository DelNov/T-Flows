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
  class(Amg_Type), target :: Amg
  integer                 :: level, ivstar
  real                    :: a(:), u(:), f(:)
  integer                 :: ia(:), ja(:)
  integer                 :: iw(:)
!-----------------------------------[locals]-----------------------------------!
  real                         :: fac, s1, s2, sa
  integer                      :: i, iaux, j, n
  real,    contiguous, pointer :: lev_a(:), lev_u(:), lev_f(:)
  integer, contiguous, pointer :: lev_ia(:), lev_ja(:)
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  if(ivstar .ne. 1) return

  !-----------------------------------------!
  !   Computation of scaling factor "fac"   !
  !-----------------------------------------!

#ifdef AMG_USE_OLD_LOOP

  ! See comment in source "Coarsening.f90" at line 180
  iaux = ia(Amg % imax(level)+1)
  ia(Amg % imax(level)+1) = iw(Amg % iminw(level))

  s1 = 0.0
  s2 = 0.0
  do i = Amg % imin(level), Amg % imax(level)
    sa = 0.0
    do j = ia(i), ia(i+1) - 1
      sa = sa+a(j)*u(ja(j))
    end do
    s1 = s1+u(i)*f(i)
    s2 = s2+u(i)*sa
  end do

  fac = 1.0
  if(s2 .ne. 0.0) fac = s1/s2

  !-------------!
  !   Scaling   !
  !-------------!
  do i = Amg % imin(level), Amg % imax(level)
    u(i) = u(i)*fac
  end do
  ia(Amg % imax(level)+1) = iaux
#endif

#ifdef AMG_USE_NEW_LOOP
  n      =  Amg % lev(level) % n
  lev_a  => Amg % lev(level) % a
  lev_u  => Amg % lev(level) % u
  lev_f  => Amg % lev(level) % f
  lev_ia => Amg % lev(level) % ia
  lev_ja => Amg % lev(level) % ja

  call Amg % Update_U_And_F_At_Level(level, vec_u=u, vec_f=f)

  s1 = 0.0
  s2 = 0.0
  do i = 1, n
    sa = 0.0
    do j = lev_ia(i), lev_ia(i+1) - 1
      sa = sa + lev_a(j) * lev_u(lev_ja(j))
    end do
    s1 = s1 + lev_u(i) * lev_f(i)
    s2 = s2 + lev_u(i) * sa
  end do

  fac = 1.0
  if(s2 .ne. 0.0) fac = s1/s2

  !-------------!
  !   Scaling   !
  !-------------!
  do i = 1, n
    lev_u(i) = lev_u(i) * fac
  end do

  call Amg % Update_U_And_F_Globally(level, vec_u=u)
#endif

  end subroutine
