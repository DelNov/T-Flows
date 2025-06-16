!==============================================================================!
  subroutine Interpolate_Correction(Amg, level)
!------------------------------------------------------------------------------!
!   Interpolates correction from grid level+1 to grid level
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(Amg_Type), target :: Amg
  integer                 :: level
!-----------------------------------[locals]-----------------------------------!
  integer                      :: i, j, ic, if, n, n_c, nw
  real,    contiguous, pointer :: c_u(:)
  integer, contiguous, pointer :: c_ifg_1(:)
  real,    contiguous, pointer :: u(:)
  integer, contiguous, pointer :: ifg_2(:)
  integer, contiguous, pointer :: jw(:)
  real,    contiguous, pointer :: w(:)
  integer :: i_min, i_max
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  !--------------------------!
  !   c -> c contributions   !
  !--------------------------!
  n_c     =  Amg % lev(level+1) % n      ! number of cells on coarse grid
  c_u     => Amg % lev(level+1) % u      ! u on coarse grid
  c_ifg_1 => Amg % lev(level+1) % ifg_1  ! ifg_1 on coarse grid

  n     =  Amg % lev(level) % n
  nw    =  Amg % lev(level) % nw
  u     => Amg % lev(level) % u
  ifg_2 => Amg % lev(level) % ifg_2
  w     => Amg % lev(level) % w
  jw    => Amg % lev(level) % jw

  do ic = 1, n_c
    if = c_ifg_1(ic)
    u(if) = u(if) + c_u(ic)
  end do

  !--------------------------!
  !   c -> f contributions   !
  !--------------------------!
  do i = 1, nw
    if = ifg_2(i)
    do j = Amg % lev(level) % iw(i), Amg % lev(level) % iw(i+1)-1
      u(if) = u(if) + w(j) * c_u(jw(j))
    end do
  end do

  end subroutine
