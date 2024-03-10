!==============================================================================!
  subroutine Vec_D_Vec(Lin, n, dot, A, B)
!------------------------------------------------------------------------------!
!>  Front-end for calculation of vector vector dot product.
!------------------------------------------------------------------------------!
!   Note: Using intent clause here, was causing slower runs and crashes        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Linalg_Type)  :: Lin   !! parent class
  integer, intent(in) :: n     !! size of vectors
  real                :: dot   !! result of the dot product
  real                :: a(n)  !! operand vector
  real                :: b(n)  !! operand vector
!==============================================================================!

  call Lin % Vec_D_Vec_Acc(dot, n, a, b)

  end subroutine

