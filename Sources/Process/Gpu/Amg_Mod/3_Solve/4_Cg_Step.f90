!==============================================================================!
  subroutine Cg_Step(Amg, level, icgr, iter,  &
                     a, u, f, ia, ja,         &  ! defining system
                     iw,                      &
                     m)
!------------------------------------------------------------------------------!
!   Performs one step of preconditioned conjugate gradient
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(Amg_Type) :: Amg
  integer         :: level, icgr, iter
  real            :: a(:), u(:), f(:)
  integer         :: ia(:), ja(:)
  integer         :: iw(:)
  integer         :: m
!-----------------------------------[locals]-----------------------------------!
  real    :: alf, eps, s2
  integer :: i, ishift, nnu
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  if (icgr.eq.0) return

  nnu = Amg % imax(level)-Amg % imin(level)+1
  ishift = Amg % imax(m)+1-Amg % imin(level)

  !---------------------------------------!
  !   Compute most recent MG correction   !
  !---------------------------------------!
  do i = Amg % imin(level), Amg % imax(level)
    u(i) = u(i) - u(i+ishift)
  end do

  !-------------------!
  !   First CG step   !
  !-------------------!
  if(icgr.eq.1 .or. iter.le.1) then
    do i = Amg % imin(level), Amg % imax(level)
      f(i+ishift) = u(i)
    end do

  !------------------------------------------!
  !   Subsequent CG steps (only if icgr=2)   !
  !------------------------------------------!
  else
    alf = Amg % Cg_Alpha(level, s2,        &
                         a, u, f, ia, ja,  &
                         iw, m)
    do i = Amg % imin(level), Amg % imax(level)
      f(i+ishift) = u(i)+alf*f(i+ishift)
    end do
  end if

  eps = Amg % Cg_Epsilon(level, s2, a, u, f, ia, ja, iw, m)

  if(Amg % ierr .gt. 0) then
    return
  end if

  do i = Amg % imin(level), Amg % imax(level)
    u(i) = u(i+ishift) + eps * f(i+ishift)
  end do

  end subroutine
