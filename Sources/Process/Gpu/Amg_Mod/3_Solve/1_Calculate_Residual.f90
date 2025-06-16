!==============================================================================!
  subroutine Calculate_Residual(Amg, level, resl,  &
                                a, u, f, ia, ja,   &  ! defining system
                                iw)
!------------------------------------------------------------------------------!
!   Computes l2-norm of residual on level "level"
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(Amg_Type), target :: Amg
  integer                 :: level
  real                    :: resl
  real                    :: a(:), u(:), f(:)
  integer                 :: ia(:), ja(:)
  integer                 :: iw(:)
!-----------------------------------[locals]-----------------------------------!
  real                         :: s
  integer                      :: i, iaux, j, n, i_loc
  real,    contiguous, pointer :: lev_a(:), lev_u(:), lev_f(:)
  integer, contiguous, pointer :: lev_ia(:), lev_ja(:)
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  resl = 0.0

!old:  ! See comment in source "Coarsening.f90" at line 180
!old:  iaux = ia(Amg % imax(level)+1)
!old:  ia(Amg % imax(level)+1) = iw(Amg % iminw(level))
!old:
!old:   do i = Amg % imin(level), Amg % imax(level)
!old:     s = f(i)
!old:     do j = ia(i), ia(i+1) - 1
!old:       s = s - a(j) * u(ja(j))
!old:     end do
!old:     resl = resl + s*s
!old:   end do

  n      =  Amg % lev(level) % n
  lev_a  => Amg % lev(level) % a
  lev_u  => Amg % lev(level) % u
  lev_f  => Amg % lev(level) % f
  lev_ia => Amg % lev(level) % ia
  lev_ja => Amg % lev(level) % ja

  ! Copy vectors u and f to level's storage
  ! (This is needed during the development stage)
  call Amg % Update_U_And_F_At_Level(level, vec_u=u, vec_f=f)

   resl = 0.0
   do i = 1, n
     s = lev_f(i)
     do j = lev_ia(i), lev_ia(i+1) - 1
       s = s - lev_a(j) * lev_u(lev_ja(j))
     end do
     resl = resl + s * s
   end do

  resl = sqrt(resl)

!old:  ia(Amg % imax(level)+1) = iaux

  end subroutine
