!==============================================================================!
  pure subroutine Set_Array_Range(Math, n, minv, maxv, array)
!------------------------------------------------------------------------------!
!>  Sets values in an array of size n, uniformly distributed between the
!>  bounding values minv and maxv
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Math_Type),  intent(in)  :: Math        !! parent class
  integer,           intent(in)  :: n           !! size of the array
  real,              intent(in)  :: minv, maxv  !! bounding values
  real,              intent(out) :: array(n)    !! resulting array
!-----------------------------------[Locals]-----------------------------------!
  integer :: i
  real    :: delta
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Math)
!==============================================================================!

  array(1) = minv
  array(n) = maxv
  delta = (maxv-minv) / real(n-1.0)
  do i = 1, n
    array(i) = minv + delta * (i-1)
  end do

  end subroutine
