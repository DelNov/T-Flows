!==============================================================================!
  subroutine Cg_Step(Amg, level, iter)
!------------------------------------------------------------------------------!
!   Performs one step of preconditioned conjugate gradient
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(Amg_Type), target :: Amg
  integer                 :: level, iter
!-----------------------------------[locals]-----------------------------------!
  real    :: alf, eps, s2
  integer :: i, n
  real,    contiguous, pointer :: u(:), u_b(:), f_b(:)
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  if(Amg % cycle % cg_usage .eq. AMG_NO_CG_STEPS) return

  n   =  Amg % lev(level) % n
  u   => Amg % lev(level) % u
  u_b => Amg % lev(level) % u_b
  f_b => Amg % lev(level) % f_b

  !---------------------------------------!
  !   Compute most recent MG correction   !
  !---------------------------------------!
  do i = 1, n
    u(i) = u(i) - u_b(i)
  end do

  !-------------------!
  !   First CG step   !
  !-------------------!
  if(Amg % cycle % cg_usage .eq. AMG_ONE_CG_STEP .or. iter .le. 1) then
    do i = 1, n
      f_b(i) = u(i)
    end do

  !-------------------------!
  !   Subsequent CG steps   !
  !-------------------------!
  else
    alf = Amg % Cg_Alpha(level, s2)
    do i = 1, n
      f_b(i) = u(i) + alf * f_b(i)
    end do
  end if

  eps = Amg % Cg_Epsilon(level, s2)

  if(Amg % ierr .gt. 0) then
    return
  end if

  do i = 1, n
    u(i) = u_b(i) + eps * f_b(i)
  end do

  end subroutine
