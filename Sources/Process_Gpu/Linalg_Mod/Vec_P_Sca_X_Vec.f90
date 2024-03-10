!==============================================================================!
  subroutine Vec_P_Sca_X_Vec(Lin, n, c, a, s, b)
!------------------------------------------------------------------------------!
!>  Front-end for calculation of c = a + s * b, where s is a scalar
!------------------------------------------------------------------------------!
!   Note: Using intent clause here, was causing slower runs and crashes        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Linalg_Type)  :: Lin   !! parent class
  integer, intent(in) :: n     !! size of vectors
  real                :: c(n)  !! result vector
  real                :: a(n)  !! operand vector
  real                :: s     !! operand scalar
  real                :: b(n)  !! operand vector
!==============================================================================!

  call Lin % Vec_P_Sca_X_Vec_Acc(n, c, a, s, b)

  end subroutine

