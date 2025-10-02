!==============================================================================!
  subroutine Sys_Restore(Lin, n, fn, A, b)
!------------------------------------------------------------------------------!
!>  Front-end for restoring (de-normalizing) a linear system of equations
!------------------------------------------------------------------------------!
!   Note: Using intent clause here, was causing slower runs and crashes        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Linalg_Type)        :: Lin   !! parent class
  integer, intent(in)       :: n     !! size of vectors
  real,    intent(in)       :: fn    !! factor of normalization
  type(Sparse_Type), target :: A     !! operand matrix
  real                      :: b(n)  !! right hand side vector
!-----------------------------------[Locals]-----------------------------------!
  real, pointer :: a_val(:), d_inv(:)
  real          :: fn_inv
  integer       :: i, nz
!==============================================================================!

  ! Take aliases
  nz    =  A % nonzeros
  a_val => A % val
  d_inv => A % d_inv

  !------------------------------------------------!
  !   Work out the inverse of the scaling factor   !
  !------------------------------------------------!
  fn_inv = 1.0 / fn       ! inverse scaling factor, like 1.0 / fn

  !-----------------------------------------------------!
  !   Scale the matrix and the right hand side vector   !
  !-----------------------------------------------------!

  !$tf-acc loop begin
  do i = 1, nz
    a_val(i) = a_val(i) * fn_inv
  end do
  !$tf-acc loop end

  !$tf-acc loop begin
  do i = 1, n
    b(i)     = b(i)     * fn_inv
    d_inv(i) = d_inv(i) * fn
  end do
  !$tf-acc loop end

  end subroutine

