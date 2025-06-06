!==============================================================================!
  double precision function cg_epsilon(amg, level, s2,   &
                                       a, u, f, ia, ja,  &  ! defining system
                                       iw,               &
                                       m)
!------------------------------------------------------------------------------!
!   Called from cg_step
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
  double precision :: s1, sp, sr
  integer          :: i, iaux, ishift, j
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  ishift = amg % imax(m)+1-amg % imin(level)
  s1 = 0.0d0
  s2 = 0.0d0

  ! See comment in source "coarsening.f90" at line 180
  iaux = ia(amg % imax(level)+1)
  ia(amg % imax(level)+1) = iw(amg % iminw(level))

  do i = amg % imin(level), amg % imax(level)
    sr = f(i)
    sp = 0.0d0
    do j = ia(i), ia(i+1) - 1
      sr = sr - a(j) * u(ja(j) + ishift)
      sp = sp + a(j) * f(ja(j) + ishift)
    end do
    s1 = s1 + sr * f(i+ishift)
    s2 = s2 + sp * f(i+ishift)
  end do

  ia(amg % imax(level)+1) = iaux

  ! Error exit
  if(s2 .eq. 0.0d0) then
    write(6, '(a)')  &
      ' *** error in cgeps: cg correction not defined ***'
    amg % ierr = AMG_ERR_CG_CORRECTION_UNDEF
    return
  end if

  cg_epsilon = +s1 / s2

  end function
