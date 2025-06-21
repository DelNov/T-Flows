!==============================================================================!
  real function Cg_Alpha(Amg, level, s2)
!------------------------------------------------------------------------------!
!   Called from Cg_Step, it seems to be updating alpha from the CG method
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(Amg_Type), target :: Amg
  integer                 :: level
  real                    :: s2
!-----------------------------------[locals]-----------------------------------!
  real                         :: s1, sr
  integer                      :: i, j, ij, n
  real,    contiguous, pointer :: a(:), u(:), f_b(:)
  integer, contiguous, pointer :: ia(:), ja(:)
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  s1 = 0.0

  n   =  Amg % lev(level) % n
  a   => Amg % lev(level) % a
  ia  => Amg % lev(level) % ia
  ja  => Amg % lev(level) % ja
  u   => Amg % lev(level) % u
  f_b => Amg % lev(level) % f_b

  do i = 1, n
    sr = 0.0
    do ij = ia(i), ia(i+1) - 1
      j = ja(ij)
      sr = sr + a(ij) * u(j)
    end do
    s1 = s1 + sr * f_b(i)
  end do

  Cg_Alpha = -s1 / s2

  end function
