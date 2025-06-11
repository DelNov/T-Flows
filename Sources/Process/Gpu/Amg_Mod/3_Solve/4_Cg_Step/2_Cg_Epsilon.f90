!==============================================================================!
  real     function Cg_Epsilon(Amg, level, s2,   &
                                       a, u, f, ia, ja,  &  ! defining system
                                       iw,               &
                                       m)
!------------------------------------------------------------------------------!
!   Called from Cg_Step
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
  real    :: s1, sp, sr
  integer :: i, iaux, ishift, j
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  ishift = Amg % imax(m)+1-Amg % imin(level)
  s1 = 0.0
  s2 = 0.0

  ! See comment in source "Coarsening.f90" at line 180
  iaux = ia(Amg % imax(level)+1)
  ia(Amg % imax(level)+1) = iw(Amg % iminw(level))

  do i = Amg % imin(level), Amg % imax(level)
    sr = f(i)
    sp = 0.0
    do j = ia(i), ia(i+1) - 1
      sr = sr - a(j) * u(ja(j) + ishift)
      sp = sp + a(j) * f(ja(j) + ishift)
    end do
    s1 = s1 + sr * f(i+ishift)
    s2 = s2 + sp * f(i+ishift)
  end do

  ia(Amg % imax(level)+1) = iaux

  ! Error exit
  if(s2 .eq. 0.0) then
    write(6, '(a)')  &
      ' *** error in cgeps: cg correction not defined ***'
    Amg % ierr = AMG_ERR_CG_CORRECTION_UNDEF
    return
  end if

  Cg_Epsilon = +s1 / s2

  end function
