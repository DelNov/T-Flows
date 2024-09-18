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
!-----------------------------------[Locals]-----------------------------------!
  integer :: i
!==============================================================================!

  !$acc parallel loop independent  &
  !$acc present(  &
  !$acc   c,  &
  !$acc   a,  &
  !$acc   b   &
  !$acc )
  do i = 1, n
    c(i) = a(i) + s * b(i)
  end do
  !$acc end parallel

  end subroutine

