!==============================================================================!
  subroutine Vec_P_Vec_X_Vec_Acc(Lin, n, d, a, b, c)
!------------------------------------------------------------------------------!
!>  This subroutine computes d = a + b * c on a device, without checking if
!>  variables are present on the device.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Linalg_Type) :: Lin   !! parent class
  integer            :: n     !! vector dimensions
  real               :: d(n)  !! result vector
  real               :: a(n)  !! operand vector
  real               :: b(n)  !! operand scalar
  real               :: c(n)  !! operand vector
!-----------------------------------[Locals]-----------------------------------!
  integer :: i
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Lin)
!==============================================================================!

  !$acc kernels
  do i = 1, n
    d(i) = a(i) + b(i) * c(i)
  end do
  !$acc end kernels

  end subroutine

