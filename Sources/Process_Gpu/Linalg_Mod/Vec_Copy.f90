!==============================================================================!
  subroutine Vec_Copy(Lin, n, a, b)
!------------------------------------------------------------------------------!
!>  Front-end for calculation of vector vector dot product.
!------------------------------------------------------------------------------!
!   Note: Using intent clause here, was causing slower runs and crashes        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Linalg_Type)  :: Lin   !! parent class
  integer, intent(in) :: n     !! size of vectors
  real                :: a(n)  !! operand vector
  real                :: b(n)  !! operand vector
!==============================================================================!

  call Lin % Vec_Copy_Acc(n, a, b)

  end subroutine

