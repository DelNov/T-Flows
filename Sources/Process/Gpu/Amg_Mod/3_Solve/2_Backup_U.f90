!==============================================================================!
  subroutine Backup_U(Amg, level, icgr, u, u_b)
!------------------------------------------------------------------------------!
!   Makes a back-up of the current approx. on level "level" if icgr.ne.0.
!   (This seems to try to place backup beyond the last (coarsest) level.)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(Amg_Type), target :: Amg
  integer                 :: level, icgr
  real                    :: u(:), u_b(:)
!-----------------------------------[locals]-----------------------------------!
  integer                      :: i, n
  real,    contiguous, pointer :: lev_u(:), lev_u_b(:)
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  if(icgr .eq. 0) return

#ifdef AMG_USE_OLD_LOOP
  do i = Amg % imin(level), Amg % imax(level)
    u_b(i) = u(i)
  end do
#endif

#ifdef AMG_USE_NEW_LOOP
  call Amg % Update_U_And_F_At_Level(level, vec_u=u)
  n       =  Amg % lev(level) % n
  lev_u   => Amg % lev(level) % u
  lev_u_b => Amg % lev(level) % u_b
  do i = 1, n
    lev_u_b(i) = lev_u(i)
  end do
  call Amg % Update_U_And_F_Globally(level, vec_u_b=u_b)
#endif

  end subroutine
