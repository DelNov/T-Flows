!==============================================================================!
  subroutine Mat_X_Vec(Lin, n, c, Acon, Aval, b)
!------------------------------------------------------------------------------!
!>  Front-end for calculation of sparse-matrix vector multiplication.
!------------------------------------------------------------------------------!
!   Note: Using intent clause here, was causing slower runs and crashes        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Linalg_Type)    :: Lin   !! parent class
  integer, intent(in)   :: n     !! size of vectors
  real                  :: c(n)  !! result vector
  type(Sparse_Con_Type) :: Acon  !! operand connectivity matrix
  type(Sparse_Val_Type) :: Aval  !! operand values matrix
  real                  :: b(n)  !! operand vector
!-----------------------------------[Locals]-----------------------------------!
  integer :: nz
!==============================================================================!

  nz = Acon % nonzeros

  call Lin % Mat_X_Vec_Acc(n,           &
                           nz,          &
                           c,           &
                           Aval % val,  &
                           Acon % col,  &
                           Acon % row,  &
                           b)

  end subroutine

