!==============================================================================!
  subroutine Restrict_Residuals(Amg, level_c)
!------------------------------------------------------------------------------!
!   Restricts residuals from grid level_c-1 to grid level_c
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(Amg_Type), target :: Amg
  integer                 :: level_c
!-----------------------------------[locals]-----------------------------------!
  real                         :: d
  integer                      :: i, ic, if, j, n_c, nw
  real,    contiguous, pointer :: c_f(:)
  integer, contiguous, pointer :: c_ifg_1(:)
  real,    contiguous, pointer :: f(:), u(:), a(:)
  integer, contiguous, pointer :: ifg_2(:)
  integer, contiguous, pointer :: ia(:), ja(:)
  integer, contiguous, pointer :: jw(:)
  real,    contiguous, pointer :: w(:)
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  !---------------------------------!
  !   Transfer of c-point defects   !
  !---------------------------------!

  ! Coarse level pointers and data
  n_c         =  Amg % lev(level_c) % n
  c_ifg_1 => Amg % lev(level_c) % ifg_1  ! ifg_1 on coarse grid
  c_f     => Amg % lev(level_c) % f

  ! Fine level pointers and data
  nw        =  Amg % lev(level_c-1) % nw
  a     => Amg % lev(level_c-1) % a
  ia    => Amg % lev(level_c-1) % ia
  ja    => Amg % lev(level_c-1) % ja
  u     => Amg % lev(level_c-1) % u
  ifg_2 => Amg % lev(level_c-1) % ifg_2
  f     => Amg % lev(level_c-1) % f
  w     => Amg % lev(level_c-1) % w
  jw    => Amg % lev(level_c-1) % jw

  do ic = 1, n_c
    if = c_ifg_1(ic)                   ! fine cell
    d = f(if)                          ! fine source
    do j = ia(if), ia(if+1) - 1    ! fine matrix entries, I hope
      d = d - a(j) * u(ja(j))  ! fine unknown
    end do
    c_f(ic) = d                        ! coarse source
  end do

  !---------------------------------!
  !   Transfer of f-point defects   !
  !---------------------------------!

  do i = 1, nw
    if = ifg_2(i)
    d = f(if)                          ! fine source
    do j = ia(if), ia(if+1) - 1    ! fine matrix entries, I hope
      d = d - a(j) * u(ja(j))  ! fine unknown
    end do
    do j = Amg % lev(level_c-1) % iw(i), Amg % lev(level_c-1) % iw(i+1)-1
      c_f(jw(j)) = c_f(jw(j)) + w(j) * d
    end do
  end do

  end subroutine
