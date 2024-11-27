!==============================================================================!
  subroutine Vec_M_Vec_X_Vec(Lin, n, d, a, b, c)
!------------------------------------------------------------------------------!
!>  Front-end for calculation of d = a - b * c, where all operands are vectors
!------------------------------------------------------------------------------!
!   Note: Using intent clause here, was causing slower runs and crashes        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Linalg_Type)  :: Lin   !! parent class
  integer, intent(in) :: n     !! size of vectors
  real                :: d(n)  !! result vector
  real                :: a(n)  !! operand vector
  real                :: b(n)  !! operand scalar
  real                :: c(n)  !! operand vector
!-----------------------------------[Locals]-----------------------------------!
  integer :: i
!==============================================================================!

  !$tf-acc loop begin
  do i = 1, n
    d(i) = a(i) - b(i) * c(i)
  end do
  !$tf-acc loop end

  end subroutine

