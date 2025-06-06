!==============================================================================!
  subroutine gauss_seidel_sweep(amg,  level, irel,  &
                                a, u, f, ia, ja,    &  ! defining system
                                iw, icg)
!------------------------------------------------------------------------------!
!   Performs one (partial) gauss-seidel sweep on grik level:
!
!   irel = 1 :   partial gauss-seidel sweep (only f-points)
!        = 2 :   full gauss-seidel sweep (all points)
!        = 3 :   partial gauss-seidel sweep (only c-points)
!        = 4 :   full sweep: ff -- c -- colors (highest first)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(amg_type)  :: amg
  integer          :: level, irel
  double precision :: a(:), u(:), f(:)
  integer          :: ia(:), ja(:)
  integer          :: iw(:), icg(:)
!-----------------------------------[locals]-----------------------------------!
  double precision :: s
  integer          :: i, iaux, j
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  call amg % timer_start()

  ! See comment in source "coarsening.f90" at line 180
  iaux = ia(amg % imax(level)+1)
  ia(amg % imax(level)+1) = iw(amg % iminw(level))

  if(irel.eq.1) then

    !------------------!
    !   F-relaxation   !
    !------------------!
    do i = amg % imin(level), amg % imax(level)
      if(icg(i).le.0) then
        s = f(i)
        do j = ia(i) + 1, ia(i+1) - 1
          s = s-a(j)*u(ja(j))
        end do
        u(i) = s/a(ia(i))
      end if
    end do

  else if(irel.eq.2) then

    !------------------------!
    !   Full GS relaxation   !
    !------------------------!
    do i = amg % imin(level), amg % imax(level)
      s = f(i)
      do j = ia(i) + 1, ia(i+1) - 1
        s = s-a(j)*u(ja(j))
      end do
      u(i) = s/a(ia(i))
    end do

  else if(irel.eq.3) then

    !------------------!
    !   C-relaxation   !
    !------------------!
    do i = amg % imin(level), amg % imax(level)
      if(icg(i) .gt. 0) then
        s = f(i)
        do j = ia(i) + 1, ia(i+1) - 1
          s = s-a(j)*u(ja(j))
        end do
        u(i) = s/a(ia(i))
      end if
    end do

  else if(irel.eq.4) then

    !-------------------!
    !   FF-relaxation   !
    !-------------------!
    do i = amg % imin(level), amg % imax(level)
      if(icg(i) .eq. 0) then
        s = f(i)
        do j = ia(i) + 1, ia(i+1) - 1
          s = s-a(j)*u(ja(j))
        end do
        u(i) = s/a(ia(i))
      end if
    end do

    !------------------!
    !   C-relaxation   !
    !------------------!
    do i = amg % imin(level), amg % imax(level)
      if(icg(i) .gt. 0) then
        s = f(i)
        do j = ia(i) + 1, ia(i+1) - 1
          s = s-a(j)*u(ja(j))
        end do
        u(i) = s/a(ia(i))
      end if
    end do

    !---------------------------!
    !   More color relaxation   !
    !---------------------------!
    i = amg % nstcol(level)
    do
      if(i .ge. AMG_BIG_INTEGER) exit
      s = f(i)
      do j = ia(i) + 1, ia(i+1) - 1
        s = s-a(j)*u(ja(j))
      end do
      u(i) = s/a(ia(i))
      i = -icg(i)
    end do

  else
    print *, "Ouch, didn't see that coming: irel is not in range "
    print *, " 1 to 4.  Go back to the original version of this  "
    print *, "  function and try to figure out what went wrong.  "
    print *, "  Oh yes, this function is gauss_seidel_sweep.f.   "
    stop
  end if

  call amg % timer_stop(13)

  ia(amg % imax(level)+1) = iaux

  end subroutine
