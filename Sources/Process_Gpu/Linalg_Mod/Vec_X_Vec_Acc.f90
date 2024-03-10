!==============================================================================!
  subroutine Vec_X_Vec_Acc(Lin, n, c, a, b)
!------------------------------------------------------------------------------!
!>  This subroutine computes vector vector multiplication on a device,
!>  without checking if variables are present on the device.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Linalg_Type) :: Lin         !! parent class
  integer            :: n           !! matrix and vector dimension
  real               :: c(n)        !! result vector
  real               :: a(n)        !! operand matrix values
  real               :: b(n)        !! operand vector
!-----------------------------------[Locals]-----------------------------------!
  integer :: i
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Lin)
!==============================================================================!

  !$acc kernels
  do i = 1, n
    c(i) = a(i) * b(i)
  end do
  !$acc end kernels

  end subroutine

