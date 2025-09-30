!==============================================================================!
  subroutine Gauss_Seidel_Sweep(Amg,  level, irel)
!------------------------------------------------------------------------------!
!   Performs one (partial) gauss-seidel sweep on grik level:
!
!   irel = AMG_RELAX_F_POINTS  (1):  partial Gauss-Seidel sweep (only f-points)
!        = AMG_RELAX_FULL_GS   (2):  full Gauss-Seidel sweep (all points)
!        = AMG_RELAX_C_POINTS  (3):  partial Gauss-Seidel sweep (only c-points)
!        = AMG_RELAX_MULTICOLOR(4):  full sweep: ff--c--colors (highest first)
!------------------------------------------------------------------------------!
  implicit none
!----------------------------------[Arguments]----------------------------------!
  class(Amg_Type), target :: Amg
  integer                 :: level, irel
!-----------------------------------[Locals]-----------------------------------!
  real,    contiguous, pointer :: a(:), u(:), f(:)
  integer, contiguous, pointer :: ia(:), ja(:), icg(:)
  real                         :: s
  integer                      :: i, j, ij, n
!------------------------------------[Save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  ! Fetch some pointers
  n   =  Amg % lev(level) % n
  a   => Amg % lev(level) % a
  u   => Amg % lev(level) % u
  f   => Amg % lev(level) % f
  ia  => Amg % lev(level) % ia
  ja  => Amg % lev(level) % ja
  icg => Amg % lev(level) % icg

  !------------------!
  !                  !
  !   F-relaxation   !
  !                  !
  !------------------!
  if(irel .eq. AMG_RELAX_F_POINTS) then

    do i = 1, n
      if(icg(i) .le. 0) then
        s = f(i)
        do ij = ia(i) + 1, ia(i+1) - 1
          j = ja(ij)
          s = s - a(ij) * u(j)
        end do
        u(i) = s / a(ia(i))
      end if
    end do

  !------------------------!
  !                        !
  !   Full GS relaxation   !
  !                        !
  !------------------------!
  else if(irel .eq. AMG_RELAX_FULL_GS) then

    do i = 1, n
      s = f(i)
      do ij = ia(i) + 1, ia(i+1) - 1
        j = ja(ij)
        s = s - a(ij) * u(j)
      end do
      u(i) = s / a(ia(i))
    end do

  !------------------!
  !                  !
  !   C-relaxation   !
  !                  !
  !------------------!
  else if(irel .eq. AMG_RELAX_C_POINTS) then

    do i = 1, n
      if(icg(i) .gt. 0) then
        s = f(i)
        do ij = ia(i) + 1, ia(i+1) - 1
          j = ja(ij)
          s = s - a(ij) * u(j)
        end do
        u(i) = s / a(ia(i))
      end if
    end do

  !----------------------------!
  !                            !
  !   Multi-color relaxation   !
  !                            !
  !----------------------------!
  else if(irel.eq.AMG_RELAX_MULTICOLOR) then

    !-------------------!
    !   FF-relaxation   !
    !-------------------!
    do i = 1, n
      if(icg(i) .eq. 0) then
        s = f(i)
        do ij = ia(i) + 1, ia(i+1) - 1
          j = ja(ij)
          s = s - a(ij) * u(j)
        end do
        u(i) = s / a(ia(i))
      end if
    end do

    !------------------!
    !   C-relaxation   !
    !------------------!
    do i = 1, n
      if(icg(i) .gt. 0) then
        s = f(i)
        do ij = ia(i) + 1, ia(i+1) - 1
          j = ja(ij)
          s = s - a(ij) * u(j)
        end do
        u(i) = s / a(ia(i))
      end if
    end do

    !----------------------------!
    !   Multi-color relaxation   !
    !----------------------------!
    i = Amg % lev(level) % start_of_color
    do
      if(i .ge. AMG_BIG_INTEGER) exit
      s = f(i)
      do ij = ia(i) + 1, ia(i+1) - 1
        j = ja(ij)
        s = s - a(ij) * u(j)
      end do
      u(i) = s / a(ia(i))
      i = -icg(i)
      Assert(i .gt. 0)
    end do

  else
    print *, "Ouch, didn't see that coming: irel is not in range "
    print *, " 1 to 4.  Go back to the original version of this  "
    print *, "  function and try to figure out what went wrong.  "
    print *, "  Oh yes, this function is Gauss_Seidel_Sweep.f.   "
    stop
  end if

  end subroutine
