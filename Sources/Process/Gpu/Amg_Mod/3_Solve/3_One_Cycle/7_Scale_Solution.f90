!==============================================================================!
  subroutine Scale_Solution(Amg, level, ivstar)
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
!-----------------------------------[locals]-----------------------------------!
  real                         :: fac, s1, s2, sa
  integer                      :: i, j, ij, n
  real,    contiguous, pointer :: a(:), u(:), f(:)
  integer, contiguous, pointer :: ia(:), ja(:)
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  if(ivstar .ne. 1) return

  !-----------------------------------------!
  !   Computation of scaling factor "fac"   !
  !-----------------------------------------!

  n  =  Amg % lev(level) % n
  a  => Amg % lev(level) % a
  u  => Amg % lev(level) % u
  f  => Amg % lev(level) % f
  ia => Amg % lev(level) % ia
  ja => Amg % lev(level) % ja

  s1 = 0.0
  s2 = 0.0
  do i = 1, n
    sa = 0.0
    do ij = ia(i), ia(i+1) - 1
      j = ja(ij)
      sa = sa + a(ij) * u(j)
    end do
    s1 = s1 + u(i) * f(i)
    s2 = s2 + u(i) * sa
  end do

  fac = 1.0
  if(s2 .ne. 0.0) fac = s1/s2

  !-------------!
  !   Scaling   !
  !-------------!
  do i = 1, n
    u(i) = u(i) * fac
  end do

  end subroutine
