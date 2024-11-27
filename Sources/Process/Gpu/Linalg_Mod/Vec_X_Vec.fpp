!==============================================================================!
  subroutine Vec_X_Vec(Lin, n, c, a, b)
!------------------------------------------------------------------------------!
!>  Front-end for calculation of vector vector multiplication.
!------------------------------------------------------------------------------!
!   Note: Using intent clause here, was causing slower runs and crashes        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Linalg_Type)  :: Lin   !! parent class
  integer, intent(in) :: n     !! size of vectors
  real                :: c(n)  !! result vector
  real                :: a(n)  !! operand vector
  real                :: b(n)  !! operand vector
!-----------------------------------[Locals]-----------------------------------!
  integer :: i
!==============================================================================!

  !$tf-acc loop begin
  do i = 1, n
    c(i) = a(i) * b(i)
  end do
  !$tf-acc loop end

  end subroutine

