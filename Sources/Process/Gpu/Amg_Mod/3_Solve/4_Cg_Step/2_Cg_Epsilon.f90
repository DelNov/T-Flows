!==============================================================================!
  real function Cg_Epsilon(Amg, level, s2)
!------------------------------------------------------------------------------!
!   Called from Cg_Step
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(Amg_Type), target :: Amg
  integer                 :: level
  real                    :: s2
!-----------------------------------[locals]-----------------------------------!
  real                         :: s1, sp, sr
  integer                      :: i, j, n
  real,    contiguous, pointer :: a(:), u_b(:), f(:), f_b(:)
  integer, contiguous, pointer :: ia(:), ja(:)
!------------------------------------[save]------------------------------------!
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
    do j = ia(i), ia(i+1) - 1
      sr = sr - a(j) * u_b(ja(j))
      sp = sp + a(j) * f_b(ja(j))
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
