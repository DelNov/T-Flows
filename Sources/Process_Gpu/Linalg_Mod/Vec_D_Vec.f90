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
!-----------------------------------[Locals]-----------------------------------!
  integer :: i
!==============================================================================!

  dot = 0.0

  !$acc parallel loop independent reduction(+: dot)  &
  !$acc present(  &
  !$acc   a,  &
  !$acc   b   &
  !$acc )
  do i = 1, n
    dot = dot + a(i) * b(i)
  end do
  !$acc end parallel

  ! ... then make a global sum over all processors.
  call Global % Sum_Real(dot)

  end subroutine

