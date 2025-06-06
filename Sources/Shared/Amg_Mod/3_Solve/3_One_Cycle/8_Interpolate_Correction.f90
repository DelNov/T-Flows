!==============================================================================!
  subroutine Interpolate_Correction(Amg, level,    &
                                    a, u, ia, ja,  &  ! defining system
                                    iw, ifg)
!------------------------------------------------------------------------------!
!   Interpolates correction from grid level+1 to grid level
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(Amg_Type)  :: Amg
  integer          :: level
  double precision :: a(:), u(:)
  integer          :: ia(:), ja(:)
  integer          :: iw(:), ifg(:)
!-----------------------------------[locals]-----------------------------------!
  integer :: i, j, ic, if, iaux
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  call Amg % timer_start()

  !--------------------------!
  !   c -> c contributions   !
  !--------------------------!
  do ic = Amg % imin(level+1), Amg % imax(level+1)
    if = ifg(ic)
    u(if) = u(if) + u(ic)
  end do

  !--------------------------!
  !   c -> f contributions   !
  !--------------------------!
  iaux = iw(Amg % imaxw(level)+1)
  iw(Amg % imaxw(level)+1) = ia(Amg % imin(level+1))

  do i = Amg % iminw(level), Amg % imaxw(level)
    if = ifg(i)
    do j = iw(i), iw(i+1) - 1
      u(if) = u(if)+a(j)*u(ja(j))
    end do
  end do

  iw(Amg % imaxw(level)+1) = iaux

  call Amg % timer_stop(11)

  end subroutine
