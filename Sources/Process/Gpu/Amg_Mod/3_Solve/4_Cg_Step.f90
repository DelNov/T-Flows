!==============================================================================!
  subroutine Cg_Step(Amg, level, icgr, iter,     &
                     a, u, u_b, f, f_b, ia, ja,  &  ! defining system
                     iw)
!------------------------------------------------------------------------------!
!   Performs one step of preconditioned conjugate gradient
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(Amg_Type), target :: Amg
  integer                 :: level, icgr, iter
  real                    :: a(:), u(:), u_b(:), f(:), f_b(:)
  integer                 :: ia(:), ja(:)
  integer                 :: iw(:)
!-----------------------------------[locals]-----------------------------------!
  real    :: alf, eps, s2
  integer :: i, n
  real,    contiguous, pointer :: lev_u(:), lev_u_b(:), lev_f_b(:)
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  if (icgr.eq.0) return

  n       =  Amg % lev(level) % n
  lev_u   => Amg % lev(level) % u
  lev_u_b => Amg % lev(level) % u_b
  lev_f_b => Amg % lev(level) % f_b

  !---------------------------------------!
  !   Compute most recent MG correction   !
  !---------------------------------------!
#ifdef AMG_USE_OLD_LOOP
  do i = Amg % imin(level), Amg % imax(level)
    u(i) = u(i) - u_b(i)
  end do
#endif

#ifdef AMG_USE_NEW_LOOP
  call Amg % Update_U_And_F_At_Level(level, vec_u=u, vec_u_b=u_b)
  do i = 1, n
    lev_u(i) = lev_u(i) - lev_u_b(i)
  end do
  call Amg % Update_U_And_F_Globally(level, vec_u=u)
#endif

  !-------------------!
  !   First CG step   !
  !-------------------!
  if(icgr.eq.1 .or. iter.le.1) then
#ifdef AMG_USE_OLD_LOOP
    do i = Amg % imin(level), Amg % imax(level)
      f_b(i) = u(i)
    end do
#endif

#ifdef AMG_USE_NEW_LOOP
    call Amg % Update_U_And_F_At_Level(level, vec_u=u)
    do i = 1, n
      lev_f_b(i) = lev_u(i)
    end do
    call Amg % Update_U_And_F_Globally(level, vec_f_b=f_b)
#endif

  !------------------------------------------!
  !   Subsequent CG steps (only if icgr=2)   !
  !------------------------------------------!
  else
    alf = Amg % Cg_Alpha(level, s2,             &
                         a, u, f, f_b, ia, ja,  &
                         iw)
#ifdef AMG_USE_OLD_LOOP
    do i = Amg % imin(level), Amg % imax(level)
      f_b(i) = u(i) + alf * f_b(i)
    end do
#endif

#ifdef AMG_USE_NEW_LOOP
    call Amg % Update_U_And_F_At_Level(level, vec_u=u, vec_f_b=f_b)
    do i = 1, n
      lev_f_b(i) = lev_u(i) + alf * lev_f_b(i)
    end do
    call Amg % Update_U_And_F_Globally(level, vec_f_b=f_b)
#endif
  end if

  eps = Amg % Cg_Epsilon(level, s2, a, u, u_b, f, f_b, ia, ja, iw)

  if(Amg % ierr .gt. 0) then
    return
  end if

#ifdef AMG_USE_OLD_LOOP
  do i = Amg % imin(level), Amg % imax(level)
    u(i) = u_b(i) + eps * f_b(i)
  end do
#endif

#ifdef AMG_USE_NEW_LOOP
  call Amg % Update_U_And_F_At_Level(level, vec_u_b=u_b, vec_f_b=f_b)
  do i = 1, n
    lev_u(i) = lev_u_b(i) + eps * lev_f_b(i)
  end do
  call Amg % Update_U_And_F_Globally(level, vec_u=u)
#endif

  end subroutine
