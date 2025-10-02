!==============================================================================!
  subroutine Calculate_Residual(Amg, level, resl)
!------------------------------------------------------------------------------!
!   Computes l2-norm of residual on level "level"
!------------------------------------------------------------------------------!
  implicit none
!----------------------------------[Arguments]----------------------------------!
  class(Amg_Type), target :: Amg
  integer                 :: level
  real                    :: resl
!-----------------------------------[Locals]-----------------------------------!
  real                         :: s
  integer                      :: i, j, ij, n
  real,    contiguous, pointer :: a(:), u(:), f(:)
  integer, contiguous, pointer :: ia(:), ja(:)
!------------------------------------[Save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  resl = 0.0

  n  =  Amg % lev(level) % n
  a  => Amg % lev(level) % a
  u  => Amg % lev(level) % u
  f  => Amg % lev(level) % f
  ia => Amg % lev(level) % ia
  ja => Amg % lev(level) % ja

  do i = 1, n
    s = f(i)
    do ij = ia(i), ia(i+1) - 1
      j = ja(ij)
      s = s - a(ij) * u(j)
    end do
    resl = resl + s * s
  end do

  resl = sqrt(resl)

  end subroutine
