!==============================================================================!
  subroutine restrict_residuals(amg, level_c,     &
                                a, u, f, ia, ja,  &
                                iw, ifg)
!------------------------------------------------------------------------------!
!   Restricts residuals from grid level_c-1 to grid level_c
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(amg_type)  :: amg
  integer          :: level_c
  double precision :: a(:), u(:), f(:)
  integer          :: ia(:), ja(:)
  integer          :: iw(:), ifg(:)
!-----------------------------------[locals]-----------------------------------!
  double precision :: d
  integer          :: i, iaux, iaux1, ic, if, j
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  call amg % timer_start()

  !---------------------------------!
  !   Transfer of c-point defects   !
  !---------------------------------!
  iaux = ia(amg % imax(level_c-1)+1)
  ia(amg % imax(level_c-1)+1) = iw(amg % iminw(level_c-1))
  iaux1 = iw(amg % imaxw(level_c-1)+1)
  iw(amg % imaxw(level_c-1)+1) = iaux

  do ic = amg % imin(level_c), amg % imax(level_c)
    if = ifg(ic)
    d = f(if)
    do j = ia(if), ia(if+1) - 1
      d = d - a(j) * u(ja(j))
    end do
    f(ic) = d
  end do

  !---------------------------------!
  !   Transfer of f-point defects   !
  !---------------------------------!
  do i = amg % iminw(level_c-1), amg % imaxw(level_c-1)
    if = ifg(i)
    d = f(if)
    do j = ia(if), ia(if+1) - 1
      d = d-a(j)*u(ja(j))
    end do
    do j = iw(i), iw(i+1)-1
      f(ja(j)) = f(ja(j))+a(j)*d
    end do
  end do
  ia(amg % imax(level_c-1)+1) = iaux
  iw(amg % imaxw(level_c-1)+1) = iaux1

  call amg % timer_stop(12)

  end subroutine
