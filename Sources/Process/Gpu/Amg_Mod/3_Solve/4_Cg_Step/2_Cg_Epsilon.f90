!==============================================================================!
  real function Cg_Epsilon(Amg, level, s2,             &
                           a, u, u_b, f, f_b, ia, ja,  &  ! define system
                           iw)
!------------------------------------------------------------------------------!
!   Called from Cg_Step
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(Amg_Type), target :: Amg
  integer                 :: level
  real                    :: s2
  real                    :: a(:), u(:), u_b(:), f(:), f_b(:)
  integer                 :: ia(:), ja(:)
  integer                 :: iw(:)
!-----------------------------------[locals]-----------------------------------!
  real                         :: s1, sp, sr
  integer                      :: i, iaux, j, n
  real,    contiguous, pointer :: lev_a(:), lev_u_b(:), lev_f(:), lev_f_b(:)
  integer, contiguous, pointer :: lev_ia(:), lev_ja(:)
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  s1 = 0.0
  s2 = 0.0

  ! See comment in source "Coarsening.f90" at line 180
#ifdef AMG_USE_OLD_LOOP
  iaux = ia(Amg % imax(level)+1)
  ia(Amg % imax(level)+1) = iw(Amg % iminw(level))

  do i = Amg % imin(level), Amg % imax(level)
    sr = f(i)
    sp = 0.0
    do j = ia(i), ia(i+1) - 1
      sr = sr - a(j) * u_b(ja(j))
      sp = sp + a(j) * f_b(ja(j))
    end do
    s1 = s1 + sr * f_b(i)
    s2 = s2 + sp * f_b(i)
  end do

  ia(Amg % imax(level)+1) = iaux
#endif

#ifdef AMG_USE_NEW_LOOP
  call Amg % Update_U_And_F_At_Level(level, vec_u_b=u_b, vec_f=f, vec_f_b=f_b)
  n       =  Amg % lev(level) % n
  lev_a   => Amg % lev(level) % a
  lev_ia  => Amg % lev(level) % ia
  lev_ja  => Amg % lev(level) % ja
  lev_u_b => Amg % lev(level) % u_b
  lev_f   => Amg % lev(level) % f
  lev_f_b => Amg % lev(level) % f_b

  do i = 1, n
    sr = lev_f(i)
    sp = 0.0
    do j = lev_ia(i), lev_ia(i+1) - 1
      sr = sr - lev_a(j) * lev_u_b(lev_ja(j))
      sp = sp + lev_a(j) * lev_f_b(lev_ja(j))
    end do
    s1 = s1 + sr * lev_f_b(i)
    s2 = s2 + sp * lev_f_b(i)
  end do
#endif

  ! Error exit
  if(s2 .eq. 0.0) then
    write(6, '(a)')  &
      ' *** error in cgeps: cg correction not defined ***'
    Amg % ierr = AMG_ERR_CG_CORRECTION_UNDEF
    return
  end if

  Cg_Epsilon = +s1 / s2

  end function
