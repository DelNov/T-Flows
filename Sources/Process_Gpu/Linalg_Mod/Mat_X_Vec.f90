!==============================================================================!
  subroutine Mat_X_Vec(Lin, n, c, A, b)
!------------------------------------------------------------------------------!
!>  Front-end for calculation of sparse-matrix vector multiplication.
!------------------------------------------------------------------------------!
!   Note: Using intent clause here, was causing slower runs and crashes        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Linalg_Type)  :: Lin   !! parent class
  integer, intent(in) :: n     !! size of vectors
  real                :: c(n)  !! result vector
  type(Sparse_Type)   :: A     !! operand matrix
  real                :: b(n)  !! operand vector
!-----------------------------------[Locals]-----------------------------------!
  integer :: nz
!==============================================================================!

  nz = A % nonzeros

  call Lin % Mat_X_Vec_Acc(n,        &
                           nz,       &
                           c,        &
                           A % val,  &
                           A % col,  &
                           A % row,  &
                           b)

  end subroutine

