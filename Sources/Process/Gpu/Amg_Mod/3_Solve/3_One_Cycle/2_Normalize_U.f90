!==============================================================================!
  subroutine Normalize_U(Amg, level, u)
!------------------------------------------------------------------------------!
!   Normalizes u on level "level" if rowsum=0 (last component =0)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(Amg_Type), target :: Amg
  integer                 :: level
  real                    :: u(:)
!-----------------------------------[locals]-----------------------------------!
  real                      :: fac
  integer                   :: i, n
  real, contiguous, pointer :: lev_u(:)
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  if(Amg % irow0 .eq. AMG_NON_SINGULAR_MATRIX) return

  call Amg % Update_U_And_F_At_Level(level, vec_u=u)

  n      =  Amg % lev(level) % n
  lev_u  => Amg % lev(level) % u

  fac = lev_u(n)

  do i = 1, n
    lev_u(i) = lev_u(i) - fac
  end do

  call Amg % Update_U_And_F_Globally(level, vec_u=u)

  end subroutine
