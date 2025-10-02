!==============================================================================!
  subroutine Update_U_And_F_At_Level(Amg, level, vec_u, vec_u_b,  &
                                                 vec_f, vec_f_b,  &
                                                 for_real)
!------------------------------------------------------------------------------!
  implicit none
!----------------------------------[Arguments]----------------------------------!
  class(Amg_Type) :: Amg
  integer         :: level
  real,  optional :: vec_u(:)
  real,  optional :: vec_u_b(:)
  real,  optional :: vec_f(:)
  real,  optional :: vec_f_b(:)
  logical, optional :: for_real
!-----------------------------------[Locals]-----------------------------------!
  integer :: i, i_loc
!------------------------------------[Save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  if(.not. present(for_real)) return

  ! Copy vector u level's storage
  if(present(vec_u)) then
    do i = Amg % imin(level), Amg % imax(level)
      i_loc = i - Amg % imin(level) + 1
      Amg % lev(level) % u(i_loc) = vec_u(i)
    end do
  end if

  ! Copy vector u_b level's storage
  if(present(vec_u_b)) then
    do i = Amg % imin(level), Amg % imax(level)
      i_loc = i - Amg % imin(level) + 1
      Amg % lev(level) % u_b(i_loc) = vec_u_b(i)
    end do
  end if

  ! Copy vectors f to level's storage
  if(present(vec_f)) then
    do i = Amg % imin(level), Amg % imax(level)
      i_loc = i - Amg % imin(level) + 1
      Amg % lev(level) % f(i_loc) = vec_f(i)
    end do
  end if

  ! Copy vectors f_b to level's storage
  if(present(vec_f_b)) then
    do i = Amg % imin(level), Amg % imax(level)
      i_loc = i - Amg % imin(level) + 1
      Amg % lev(level) % f_b(i_loc) = vec_f_b(i)
    end do
  end if

  end subroutine
