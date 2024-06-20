!==============================================================================!
  subroutine Vec_P_Sca_X_Vec_Acc(Lin, n, c, a, s, b)
!------------------------------------------------------------------------------!
!>  This subroutine computes c = a + s * b on a device, without checking if
!>  variables are present on the device.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Linalg_Type) :: Lin   !! parent class
  integer            :: n     !! vector dimensions
  real               :: c(n)  !! result vector
  real               :: a(n)  !! operand vector
  real               :: s     !! operand scalar
  real               :: b(n)  !! operand vector
!-----------------------------------[Locals]-----------------------------------!
  integer :: i
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Lin)
!==============================================================================!

  !$acc kernels present(a, b, c)
  do i = 1, n
    c(i) = a(i) + s * b(i)
  end do
  !$acc end kernels

  end subroutine

