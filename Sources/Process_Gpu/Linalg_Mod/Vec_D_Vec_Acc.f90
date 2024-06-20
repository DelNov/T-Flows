!==============================================================================!
  subroutine Vec_D_Vec_Acc(Lin, dot, n, a, b)
!------------------------------------------------------------------------------!
!>  This subroutine computes vector vector dot product on a device,
!>  without checking if variables are present on the device.
!------------------------------------------------------------------------------!
!   Notes:                                                                     !
!                                                                              !
!   * Using intent clause here, was causing slower runs.                       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Linalg_Type) :: Lin   !! parent class
  real               :: dot   !! result of the dot product
  integer            :: n     !! vector dimensions
  real               :: a(n)  !! operand vector
  real               :: b(n)  !! operand vector
!-----------------------------------[Locals]-----------------------------------!
  integer :: i
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Lin)
!==============================================================================!

  !$acc kernels present(a, b)
  dot = 0
  do i = 1, n
    dot = dot + a(i) * b(i)
  end do
  !$acc end kernels

  end subroutine

