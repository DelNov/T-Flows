!==============================================================================!
  subroutine Set_Array_Range(Math, n, minv, maxv, array)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Math_Type)     :: Math
  integer, intent(in)  :: n
  real   , intent(in)  :: minv, maxv
  real,    intent(out) :: array(n)
!-----------------------------------[Locals]-----------------------------------!
  integer :: i
  real    :: delta
!==============================================================================!

  array(1) = minv
  array(n) = maxv
  delta = (maxv-minv) / real(n-1.0)
  do i = 1, n
    array(i) = minv + delta * (i-1)
  end do

  end subroutine
