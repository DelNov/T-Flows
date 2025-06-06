!==============================================================================!
  subroutine interpolate_correction(amg, level,    &
                                    a, u, ia, ja,  &  ! defining system
                                    iw, ifg)
!------------------------------------------------------------------------------!
!   Interpolates correction from grid level+1 to grid level
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(amg_type)  :: amg
  integer          :: level
  double precision :: a(:), u(:)
  integer          :: ia(:), ja(:)
  integer          :: iw(:), ifg(:)
!-----------------------------------[locals]-----------------------------------!
  integer :: i, j, ic, if, iaux
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  call amg % timer_start()

  !--------------------------!
  !   c -> c contributions   !
  !--------------------------!
  do ic = amg % imin(level+1), amg % imax(level+1)
    if = ifg(ic)
    u(if) = u(if) + u(ic)
  end do

  !--------------------------!
  !   c -> f contributions   !
  !--------------------------!
  iaux = iw(amg % imaxw(level)+1)
  iw(amg % imaxw(level)+1) = ia(amg % imin(level+1))

  do i = amg % iminw(level), amg % imaxw(level)
    if = ifg(i)
    do j = iw(i), iw(i+1) - 1
      u(if) = u(if)+a(j)*u(ja(j))
    end do
  end do

  iw(amg % imaxw(level)+1) = iaux

  call amg % timer_stop(11)

  end subroutine
