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
!-----------------------------------[Locals]-----------------------------------!
  integer :: i
!==============================================================================!

  !$tf-acc loop begin
  do i = 1, n
    a(i) = b(i)
  end do
  !$tf-acc loop end

  end subroutine

