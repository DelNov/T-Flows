!==============================================================================!
  real function Cg_Alpha(Amg, level, s2,        &
                         a, u, f, f_b, ia, ja,  &  ! defining system
                         iw)
!------------------------------------------------------------------------------!
!   Called from Cg_Step, it seems to be updating alpha from the CG method
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(Amg_Type), target :: Amg
  integer                 :: level
  real                    :: s2
  real                    :: a(:), u(:), f(:), f_b(:)
  integer                 :: ia(:), ja(:)
  integer                 :: iw(:)
!-----------------------------------[locals]-----------------------------------!
  real                         :: s1, sr
  integer                      :: i, iaux, j, n
  real,    contiguous, pointer :: lev_a(:), lev_u(:), lev_f_b(:)
  integer, contiguous, pointer :: lev_ia(:), lev_ja(:)
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  s1 = 0.0

  ! See comment in source "Coarsening.f90" at line 180
#ifdef AMG_USE_OLD_LOOP
  iaux = ia(Amg % imax(level)+1)
  ia(Amg % imax(level)+1) = iw(Amg % iminw(level))

  do i = Amg % imin(level), Amg % imax(level)
    sr = 0.0
    do j = ia(i), ia(i+1) - 1
      sr = sr + a(j) * u(ja(j))
    end do
    s1 = s1 + sr * f_b(i)
  end do

  Cg_Alpha = -s1 / s2

  ia(Amg % imax(level)+1) = iaux
#endif

#ifdef AMG_USE_NEW_LOOP
  call Amg % Update_U_And_F_At_Level(level, vec_u=u, vec_f_b=f_b)
  n       =  Amg % lev(level) % n
  lev_a   => Amg % lev(level) % a
  lev_ia  => Amg % lev(level) % ia
  lev_ja  => Amg % lev(level) % ja
  lev_u   => Amg % lev(level) % u
  lev_f_b => Amg % lev(level) % f_b
  do i = 1, n
    sr = 0.0
    do j = lev_ia(i), lev_ia(i+1) - 1
      sr = sr + lev_a(j) * lev_u(lev_ja(j))
    end do
    s1 = s1 + sr * lev_f_b(i)
  end do

  Cg_Alpha = -s1 / s2

#endif

  end function
