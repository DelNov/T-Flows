!==============================================================================!
  real function Cg_Epsilon(Amg, level, s2)
!------------------------------------------------------------------------------!
!   Called from Cg_Step
!------------------------------------------------------------------------------!
  implicit none
!----------------------------------[Arguments]----------------------------------!
  class(Amg_Type), target :: Amg
  integer                 :: level
  real                    :: s2
!-----------------------------------[Locals]-----------------------------------!
  real                         :: s1, sp, sr
  integer                      :: i, j, ij, n
  real,    contiguous, pointer :: a(:), u_b(:), f(:), f_b(:)
  integer, contiguous, pointer :: ia(:), ja(:)
!------------------------------------[Save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  s1 = 0.0
  s2 = 0.0

  n   =  Amg % lev(level) % n
  a   => Amg % lev(level) % a
  ia  => Amg % lev(level) % ia
  ja  => Amg % lev(level) % ja
  u_b => Amg % lev(level) % u_b
  f   => Amg % lev(level) % f
  f_b => Amg % lev(level) % f_b

  do i = 1, n
    sr = f(i)
    sp = 0.0
    do ij = ia(i), ia(i+1) - 1
      j = ja(ij)
      sr = sr - a(ij) * u_b(j)
      sp = sp + a(ij) * f_b(j)
    end do
    s1 = s1 + sr * f_b(i)
    s2 = s2 + sp * f_b(i)
  end do

  ! Error exit
  if(s2 .eq. 0.0) then
    write(6, '(a)')  &
      ' *** error in cgeps: cg correction not defined ***'
    Amg % ierr = AMG_ERR_CG_CORRECTION_UNDEF
    return
  end if

  Cg_Epsilon = +s1 / s2

  end function
