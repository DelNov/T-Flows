!==============================================================================!
  subroutine Restrict_Residuals(Amg, level_c,     &
                                a, u, f, ia, ja,  &
                                iw, ifg)
!------------------------------------------------------------------------------!
!   Restricts residuals from grid level_c-1 to grid level_c
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(Amg_Type), target :: Amg
  integer                 :: level_c
  real                    :: a(:), u(:), f(:)
  integer                 :: ia(:), ja(:)
  integer                 :: iw(:), ifg(:)
!-----------------------------------[locals]-----------------------------------!
  real                         :: d
  integer                      :: i, iaux, iaux1, ic, if, j, n_c, nw
  real,    contiguous, pointer :: lev_c_f(:)
  integer, contiguous, pointer :: lev_c_ifg_1(:)
  real,    contiguous, pointer :: lev_f(:), lev_u(:), lev_a(:)
  integer, contiguous, pointer :: lev_ifg_2(:)
  integer, contiguous, pointer :: lev_ia(:), lev_ja(:)
  integer, contiguous, pointer :: lev_jw(:)
  real,    contiguous, pointer :: lev_w(:)
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  !---------------------------------!
  !   Transfer of c-point defects   !
  !---------------------------------!
#ifdef AMG_USE_OLD_LOOP
  iaux = ia(Amg % imax(level_c-1)+1)
  ia(Amg % imax(level_c-1)+1) = iw(Amg % iminw(level_c-1))
  iaux1 = iw(Amg % imaxw(level_c-1)+1)
  iw(Amg % imaxw(level_c-1)+1) = iaux

  do ic = Amg % imin(level_c), Amg % imax(level_c)
    if = ifg(ic)
    Assert(level_c == Amg % Level_Of_Cell(if)+1)
    d = f(if)
    do j = ia(if), ia(if+1) - 1
      d = d - a(j) * u(ja(j))
    end do
    f(ic) = d
  end do
#endif

  ! Coarse level pointers and data
  n_c         =  Amg % lev(level_c) % n
  lev_c_ifg_1 => Amg % lev(level_c) % ifg_1  ! ifg_1 on coarse grid
  lev_c_f     => Amg % lev(level_c) % f

  ! Fine level pointers and data
  nw        =  Amg % lev(level_c-1) % nw
  lev_a     => Amg % lev(level_c-1) % a
  lev_ia    => Amg % lev(level_c-1) % ia
  lev_ja    => Amg % lev(level_c-1) % ja
  lev_u     => Amg % lev(level_c-1) % u
  lev_ifg_2 => Amg % lev(level_c-1) % ifg_2
  lev_f     => Amg % lev(level_c-1) % f
  lev_w     => Amg % lev(level_c-1) % w
  lev_jw    => Amg % lev(level_c-1) % jw

#ifdef AMG_USE_NEW_LOOP
  call Amg % Update_U_And_F_At_Level(level_c-1, vec_u=u)  ! fine level

  do ic = 1, n_c
    if = lev_c_ifg_1(ic)                   ! fine cell
    d = lev_f(if)                          ! fine source
    do j = lev_ia(if), lev_ia(if+1) - 1    ! fine matrix entries, I hope
      d = d - lev_a(j) * lev_u(lev_ja(j))  ! fine unknown
    end do
    lev_c_f(ic) = d                        ! coarse source
  end do

  call Amg % Update_U_And_F_Globally(level_c, vec_f=f)  ! coarse level
#endif

  !---------------------------------!
  !   Transfer of f-point defects   !
  !---------------------------------!
#ifdef AMG_USE_OLD_LOOP
  do i = Amg % iminw(level_c-1), Amg % imaxw(level_c-1)
    if = ifg(i)
    Assert(level_c-1 == Amg % Level_Of_Cell(if))
    d = f(if)
    do j = ia(if), ia(if+1) - 1
      d = d - a(j) * u(ja(j))
    end do
    do j = iw(i), iw(i+1)-1
      Assert(level_c == Amg % Level_Of_Cell(ja(j)))
      f(ja(j)) = f(ja(j)) + a(j) * d
    end do
  end do

  ia(Amg % imax(level_c-1)+1) = iaux
  iw(Amg % imaxw(level_c-1)+1) = iaux1
#endif

#ifdef AMG_USE_NEW_LOOP
  call Amg % Update_U_And_F_At_Level(level_c-1, vec_u=u, vec_f=f)  ! fine level

  do i = 1, nw
    if = lev_ifg_2(i)
    d = lev_f(if)                          ! fine source
    do j = lev_ia(if), lev_ia(if+1) - 1    ! fine matrix entries, I hope
      d = d - lev_a(j) * lev_u(lev_ja(j))  ! fine unknown
    end do
    do j = Amg % lev(level_c-1) % iw(i), Amg % lev(level_c-1) % iw(i+1)-1
      lev_c_f(lev_jw(j)) = lev_c_f(lev_jw(j)) + lev_w(j) * d
    end do
  end do

  call Amg % Update_U_And_F_Globally(level_c, vec_f=f)  ! coarse level
#endif

  end subroutine
