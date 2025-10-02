!==============================================================================!
  subroutine Interpolate_Correction(Amg, level)
!------------------------------------------------------------------------------!
!   Transfers corrections from a coarse grid (level+1) back to a fine grid (level).
!------------------------------------------------------------------------------!
  implicit none
!----------------------------------[Arguments]----------------------------------!
  class(Amg_Type), target :: Amg
  integer                 :: level
!-----------------------------------[Locals]-----------------------------------!
  integer                      :: i, k, ic, if, n, n_c, nw
  real,    contiguous, pointer :: c_u(:)
  integer, contiguous, pointer :: fine_index_direct(:)
  real,    contiguous, pointer :: u(:)
  integer, contiguous, pointer :: fine_index_weighted(:)
  integer, contiguous, pointer :: coarse_index_weighted(:)
  real,    contiguous, pointer :: w(:)
  integer                      :: i_min, i_max
!------------------------------------[Save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  !--------------------------!
  !   c -> c contributions   !
  !--------------------------!
  n_c               =  Amg % lev(level+1) % n      ! unknows on coarse grid
  c_u               => Amg % lev(level+1) % u      ! u on coarse grid
  fine_index_direct => Amg % lev(level+1) % fine_index_direct

  n                     =  Amg % lev(level) % n
  nw                    =  Amg % lev(level) % nw
  u                     => Amg % lev(level) % u
  fine_index_weighted   => Amg % lev(level) % fine_index_weighted
  w                     => Amg % lev(level) % w
  coarse_index_weighted => Amg % lev(level) % coarse_index_weighted

  do ic = 1, n_c
    if = fine_index_direct(ic)
    u(if) = u(if) + c_u(ic)      ! direct injection
  end do

  !--------------------------!
  !   c -> f contributions   !
  !--------------------------!
  do i = 1, nw
    if = fine_index_weighted(i)
    do k = Amg % lev(level) % iw(i), Amg % lev(level) % iw(i+1)-1
      u(if) = u(if) + w(k) * c_u(coarse_index_weighted(k))  ! weighted transfer
    end do
  end do

  end subroutine
