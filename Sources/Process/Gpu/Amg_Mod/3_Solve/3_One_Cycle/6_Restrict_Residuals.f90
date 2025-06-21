!==============================================================================!
  subroutine Restrict_Residuals(Amg, level_c, level)
!------------------------------------------------------------------------------!
!   Transfers residuals from a fine grid (level) to a coarse grid (level+1).
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(Amg_Type), target :: Amg
  integer                 :: level_c, level  ! level is level_c - 1
!-----------------------------------[locals]-----------------------------------!
  real                         :: d
  integer                      :: i, ic, if, j, n_c, nw
  real,    contiguous, pointer :: c_f(:)
  integer, contiguous, pointer :: fine_index_direct(:)
  real,    contiguous, pointer :: f(:), u(:), a(:)
  integer, contiguous, pointer :: fine_index_weighted(:)
  integer, contiguous, pointer :: ia(:), ja(:)
  integer, contiguous, pointer :: coarse_index_weighted(:)
  real,    contiguous, pointer :: w(:)
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  !---------------------------------!
  !   Transfer of c-point defects   !
  !---------------------------------!

  ! Coarse level pointers and data
  n_c               =  Amg % lev(level+1) % n
  fine_index_direct => Amg % lev(level+1) % fine_index_direct
  c_f               => Amg % lev(level+1) % f

  ! Fine level pointers and data
  nw                    =  Amg % lev(level) % nw
  a                     => Amg % lev(level) % a
  ia                    => Amg % lev(level) % ia
  ja                    => Amg % lev(level) % ja
  u                     => Amg % lev(level) % u
  fine_index_weighted   => Amg % lev(level) % fine_index_weighted
  f                     => Amg % lev(level) % f
  w                     => Amg % lev(level) % w
  coarse_index_weighted => Amg % lev(level) % coarse_index_weighted

  do ic = 1, n_c
    if = fine_index_direct(ic)    ! fine cell
    d = f(if)                     ! fine source
    do j = ia(if), ia(if+1) - 1   ! fine matrix entries, I hope
      d = d - a(j) * u(ja(j))     ! fine unknown
    end do
    c_f(ic) = d                   ! direct injection
  end do

  !---------------------------------!
  !   Transfer of f-point defects   !
  !---------------------------------!

  do i = 1, nw
    if = fine_index_weighted(i)
    d = f(if)                             ! fine source
    do j = ia(if), ia(if+1) - 1           ! fine matrix entries, I hope
      d = d - a(j) * u(ja(j))             ! fine unknown
    end do
    do j = Amg % lev(level) % iw(i), Amg % lev(level) % iw(i+1)-1
      c_f(coarse_index_weighted(j))  &
        = c_f(coarse_index_weighted(j)) + w(j) * d  ! weighted transfer
    end do
  end do

  end subroutine
