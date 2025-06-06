!==============================================================================!
  subroutine cg_step(amg, level, icgr, iter,  &
                     a, u, f, ia, ja,         &  ! defining system
                     iw,                      &
                     m)
!------------------------------------------------------------------------------!
!   Performs one step of preconditioned conjugate gradient
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(amg_type)  :: amg
  integer          :: level, icgr, iter
  double precision :: a(:), u(:), f(:)
  integer          :: ia(:), ja(:)
  integer          :: iw(:)
  integer          :: m
!-----------------------------------[locals]-----------------------------------!
  double precision :: alf, eps, s2
  integer          :: i, ishift, nnu
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  if (icgr.eq.0) return

  call amg % timer_start()

  nnu = amg % imax(level)-amg % imin(level)+1
  ishift = amg % imax(m)+1-amg % imin(level)

  !---------------------------------------!
  !   Compute most recent MG correction   !
  !---------------------------------------!
  do i = amg % imin(level), amg % imax(level)
    u(i) = u(i) - u(i+ishift)
  end do

  !-------------------!
  !   First CG step   !
  !-------------------!
  if(icgr.eq.1 .or. iter.le.1) then
    do i = amg % imin(level), amg % imax(level)
      f(i+ishift) = u(i)
    end do

  !------------------------------------------!
  !   Subsequent CG steps (only if icgr=2)   !
  !------------------------------------------!
  else
    alf = amg % cg_alpha(level, s2,        &
                         a, u, f, ia, ja,  &
                         iw, m)
    do i = amg % imin(level), amg % imax(level)
      f(i+ishift) = u(i)+alf*f(i+ishift)
    end do
  endif

  eps = amg % cg_epsilon(level, s2, a, u, f, ia, ja, iw, m)

  if(amg % ierr .gt. 0) then
    call amg % timer_stop(16)
    return
  endif

  do i = amg % imin(level), amg % imax(level)
    u(i) = u(i+ishift) + eps * f(i+ishift)
  end do

  call amg % timer_stop(16)

  end subroutine
