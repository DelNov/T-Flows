!==============================================================================!
  subroutine Update_U_And_F_At_Level(Amg, level, vec_u, vec_f)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(Amg_Type) :: Amg
  integer         :: level
  real,  optional :: vec_u(:)
  real,  optional :: vec_f(:)
!-----------------------------------[locals]-----------------------------------!
  integer :: i, i_loc
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  ! Copy vector u level's storage
  if(present(vec_u)) then
    do i = Amg % imin(level), Amg % imax(level)
      i_loc = i - Amg % imin(level) + 1
      Amg % lev(level) % u(i_loc) = vec_u(i)
    end do
  end if

  ! Copy vectors f to level's storage
  if(present(vec_f)) then
    do i = Amg % imin(level), Amg % imax(level)
      i_loc = i - Amg % imin(level) + 1
      Amg % lev(level) % f(i_loc) = vec_f(i)
    end do
  end if

  end subroutine
