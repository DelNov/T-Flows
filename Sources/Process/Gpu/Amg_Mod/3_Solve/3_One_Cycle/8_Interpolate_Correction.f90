!==============================================================================!
  subroutine Interpolate_Correction(Amg, level,    &
                                    a, u, ia, ja,  &  ! defining system
                                    iw, ifg)
!------------------------------------------------------------------------------!
!   Interpolates correction from grid level+1 to grid level
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(Amg_Type), target :: Amg
  integer                 :: level
  real                    :: a(:), u(:)
  integer                 :: ia(:), ja(:)
  integer                 :: iw(:), ifg(:)
!-----------------------------------[locals]-----------------------------------!
  integer                      :: i, j, ic, if, n, n_c, nw
  real,    contiguous, pointer :: lev_c_u(:)
  integer, contiguous, pointer :: lev_c_ifg_1(:)
  real,    contiguous, pointer :: lev_u(:)
  integer, contiguous, pointer :: lev_ifg_2(:)
  integer, contiguous, pointer :: lev_jw(:)
  real,    contiguous, pointer :: lev_w(:)
  integer :: i_min, i_max
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  !--------------------------!
  !   c -> c contributions   !
  !--------------------------!
  n_c         =  Amg % lev(level+1) % n    ! number of cells on coarse grid
  lev_c_u     => Amg % lev(level+1) % u    ! u on coarse grid
  lev_c_ifg_1 => Amg % lev(level+1) % ifg_1  ! ifg_1 on coarse grid

  n         =  Amg % lev(level) % n
  nw        =  Amg % lev(level) % nw
  lev_u     => Amg % lev(level) % u
  lev_ifg_2 => Amg % lev(level) % ifg_2
  lev_w     => Amg % lev(level) % w
  lev_jw    => Amg % lev(level) % jw

  call Amg % Update_U_And_F_At_Level(level+1, vec_u=u)  ! coarse level

  do ic = 1, n_c
    if = lev_c_ifg_1(ic)
    lev_u(if) = lev_u(if) + lev_c_u(ic)
  end do

  call Amg % Update_U_And_F_Globally(level, vec_u=u)  ! fine level

  !--------------------------!
  !   c -> f contributions   !
  !--------------------------!
  call Amg % Update_U_And_F_At_Level(level+1, vec_u=u)  ! coarse level

  do i = 1, nw
    if = lev_ifg_2(i)
    do j = Amg % lev(level) % iw(i), Amg % lev(level) % iw(i+1)-1
      lev_u(if) = lev_u(if) + lev_w(j) * lev_c_u(lev_jw(j))
    end do
  end do

  call Amg % Update_U_And_F_Globally(level, vec_u=u)  ! fine level

  end subroutine
