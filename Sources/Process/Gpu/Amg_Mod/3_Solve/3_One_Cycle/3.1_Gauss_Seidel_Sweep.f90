!==============================================================================!
  subroutine Gauss_Seidel_Sweep(Amg,  level, irel,  &
                                a, u, f, ia, ja,    &  ! defining system
                                iw, icg)
!------------------------------------------------------------------------------!
!   Performs one (partial) gauss-seidel sweep on grik level:
!
!   irel = AMG_RELAX_F_POINTS  (1):  partial Gauss-Seidel sweep (only f-points)
!        = AMG_RELAX_FULL_GS   (2):  full Gauss-Seidel sweep (all points)
!        = AMG_RELAX_C_POINTS  (3):  partial Gauss-Seidel sweep (only c-points)
!        = AMG_RELAX_MULTICOLOR(4):  full sweep: ff--c--colors (highest first)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(Amg_Type), target :: Amg
  integer                 :: level, irel
  real                    :: a(:), u(:), f(:)
  integer                 :: ia(:), ja(:)
  integer                 :: iw(:), icg(:)
!-----------------------------------[locals]-----------------------------------!
  real,    contiguous, pointer :: lev_a(:), lev_u(:), lev_f(:)
  integer, contiguous, pointer :: lev_ia(:), lev_ja(:), lev_icg(:)
  real                         :: s
  integer                      :: i, j, n
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  ! Fetch some pointers
  n       =  Amg % lev(level) % n
  lev_a   => Amg % lev(level) % a
  lev_u   => Amg % lev(level) % u
  lev_f   => Amg % lev(level) % f
  lev_ia  => Amg % lev(level) % ia
  lev_ja  => Amg % lev(level) % ja
  lev_icg => Amg % lev(level) % icg

  !------------------!
  !                  !
  !   F-relaxation   !
  !                  !
  !------------------!
  if(irel .eq. AMG_RELAX_F_POINTS) then

    call Amg % Update_U_And_F_At_Level(level, vec_u=u, vec_f=f)
    do i = 1, n
      if(lev_icg(i) .le. 0) then
        s = lev_f(i)
        do j = lev_ia(i) + 1, lev_ia(i+1) - 1
          s = s - lev_a(j) * lev_u(lev_ja(j))
        end do
        lev_u(i) = s / lev_a(lev_ia(i))
      end if
    end do
    call Amg % Update_U_And_F_Globally(level, vec_u=u)

  !------------------------!
  !                        !
  !   Full GS relaxation   !
  !                        !
  !------------------------!
  else if(irel .eq. AMG_RELAX_FULL_GS) then

    call Amg % Update_U_And_F_At_Level(level, vec_u=u, vec_f=f)
    do i = 1, n
      s = lev_f(i)
      do j = lev_ia(i) + 1, lev_ia(i+1) - 1
        s = s - lev_a(j) * lev_u(lev_ja(j))
      end do
      lev_u(i) = s / lev_a(lev_ia(i))
    end do
    call Amg % Update_U_And_F_Globally(level, vec_u=u)

  !------------------!
  !                  !
  !   C-relaxation   !
  !                  !
  !------------------!
  else if(irel .eq. AMG_RELAX_C_POINTS) then

    call Amg % Update_U_And_F_At_Level(level, vec_u=u, vec_f=f)
    do i = 1, n
      if(lev_icg(i) .gt. 0) then
        s = lev_f(i)
        do j = lev_ia(i)+1, lev_ia(i+1)-1
          s = s - lev_a(j) * lev_u(lev_ja(j))
        end do
        lev_u(i) = s / lev_a(lev_ia(i))
      end if
    end do
    call Amg % Update_U_And_F_Globally(level, vec_u=u)

  !----------------------------!
  !                            !
  !   Multi-color relaxation   !
  !                            !
  !----------------------------!
  else if(irel.eq.AMG_RELAX_MULTICOLOR) then

    !-------------------!
    !   FF-relaxation   !
    !-------------------!
    call Amg % Update_U_And_F_At_Level(level, vec_u=u, vec_f=f)
    do i = 1, n
      if(lev_icg(i) .eq. 0) then
        s = lev_f(i)
        do j = lev_ia(i) + 1, lev_ia(i+1) - 1
          s = s - lev_a(j) * lev_u(lev_ja(j))
        end do
        lev_u(i) = s / lev_a(lev_ia(i))
      end if
    end do
    call Amg % Update_U_And_F_Globally(level, vec_u=u)

    !------------------!
    !   C-relaxation   !
    !------------------!
    call Amg % Update_U_And_F_At_Level(level, vec_u=u, vec_f=f)
    do i = 1, n
      if(lev_icg(i) .gt. 0) then
        s = lev_f(i)
        do j = lev_ia(i) + 1, lev_ia(i+1) - 1
          s = s - lev_a(j) * lev_u(lev_ja(j))
        end do
        lev_u(i) = s / lev_a(lev_ia(i))
      end if
    end do
    call Amg % Update_U_And_F_Globally(level, vec_u=u)

    !----------------------------!
    !   Multi-color relaxation   !
    !----------------------------!
    call Amg % Update_U_And_F_At_Level(level, vec_u=u, vec_f=f)
    i = Amg % lev(level) % start_of_color
    do
      if(i .ge. AMG_BIG_INTEGER) exit
      s = lev_f(i)
      do j = lev_ia(i) + 1, lev_ia(i+1) - 1
        s = s - lev_a(j) * lev_u(lev_ja(j))
      end do
      lev_u(i) = s / lev_a(lev_ia(i))
      i = -lev_icg(i)
      Assert(i .gt. 0)
    end do
    call Amg % Update_U_And_F_Globally(level, vec_u=u)

  else
    print *, "Ouch, didn't see that coming: irel is not in range "
    print *, " 1 to 4.  Go back to the original version of this  "
    print *, "  function and try to figure out what went wrong.  "
    print *, "  Oh yes, this function is Gauss_Seidel_Sweep.f.   "
    stop
  end if

  end subroutine
