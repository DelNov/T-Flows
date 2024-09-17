!==============================================================================!
  subroutine Sys_Restore(Lin, n, fn, Acon, Aval, b)
!------------------------------------------------------------------------------!
!>  Front-end for restoring (de-normalizing) a linear system of equations
!------------------------------------------------------------------------------!
!   Note: Using intent clause here, was causing slower runs and crashes        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Linalg_Type)            :: Lin   !! parent class
  integer, intent(in)           :: n     !! size of vectors
  real,    intent(in)           :: fn    !! factor of normalization
  type(Sparse_Con_Type), target :: Acon  !! operand connectivity matrix
  type(Sparse_Val_Type), target :: Aval  !! operand values matrix
  real                          :: b(n)  !! right hand side vector
!-----------------------------------[Locals]-----------------------------------!
  real, pointer :: a_val(:), d_inv(:)
  real          :: fn_inv
  integer       :: i, nz
!==============================================================================!

  ! Take aliases
  nz    =  Acon % nonzeros
  a_val => Aval % val
  d_inv => Aval % d_inv

  !------------------------------------------------!
  !   Work out the inverse of the scaling factor   !
  !------------------------------------------------!
  fn_inv = 1.0 / fn       ! inverse scaling factor, like 1.0 / fn

  !-----------------------------------------------------!
  !   Scale the matrix and the right hand side vector   !
  !-----------------------------------------------------!

  !$acc parallel loop independent &
  !$acc present(  &
  !$acc   a_val   &
  !$acc )
  do i = 1, nz
    a_val(i) = a_val(i) * fn_inv
  end do
  !$acc end parallel

  !$acc parallel loop independent &
  !$acc present(  &
  !$acc   b,  &
  !$acc   d_inv   &
  !$acc )
  do i = 1, n
    b(i)     = b(i)     * fn_inv
    d_inv(i) = d_inv(i) * fn
  end do
  !$acc end parallel

  end subroutine

