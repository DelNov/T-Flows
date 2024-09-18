!==============================================================================!
  subroutine Sys_Normalize(Lin, n, fn, Acon, Aval, b)
!------------------------------------------------------------------------------!
!>  Front-end for scaling (normalizing) a linear system of equations
!------------------------------------------------------------------------------!
!   Note: Using intent clause here, was causing slower runs and crashes        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Linalg_Type)            :: Lin   !! parent class
  integer,  intent(in)          :: n     !! size of vectors
  real,     intent(out)         :: fn    !! factor of normalization
  type(Sparse_Con_Type), target :: Acon  !! operand connectivity matrix
  type(Sparse_Val_Type), target :: Aval  !! operand values matrix
  real                          :: b(n)  !! right hand side vector
!-----------------------------------[Locals]-----------------------------------!
  real,    pointer :: a_val(:), d_inv(:)
  integer, pointer :: a_dia(:)
  real             :: sum_a, avg_a, fn_inv
  integer          :: sum_n, i, nz
!==============================================================================!

  ! Take aliases
  nz    =  Acon % nonzeros
  a_dia => Acon % dia
  a_val => Aval % val
  d_inv => Aval % d_inv

  !-------------------------------------!
  !   Sum all the diagonal entries up   !
  !-------------------------------------!

  sum_a = 0.0

  !$acc parallel loop independent reduction(+: sum_a)  &
  !$acc present(  &
  !$acc   a_val,  &
  !$acc   a_dia   &
  !$acc )
  do i = 1, n
    sum_a = sum_a + a_val(a_dia(i))
  end do
  !$acc end parallel

  sum_n = n

  ! Make the sum a global sum
  call Global % Sum_Real(sum_a)
  call Global % Sum_Int (sum_n)  ! this is stored somewhere, check

  !-------------------------------------------------!
  !   Work out the scaling factor and its inverse   !
  !-------------------------------------------------!
  avg_a  = sum_a / sum_n
  fn     = 1.0 / avg_a    ! scaling factor
  fn_inv = 1.0 / fn       ! inverse scaling factor, like 1.0 / fn

  !-----------------------------------------------------!
  !   Scale the matrix and the right hand side vector   !
  !-----------------------------------------------------!

  !$acc parallel loop independent  &
  !$acc present(  &
  !$acc   a_val   &
  !$acc )
  do i = 1, nz
    a_val(i) = a_val(i) * fn
  end do
  !$acc end parallel

  !$acc parallel loop independent  &
  !$acc present(  &
  !$acc   b,  &
  !$acc   d_inv   &
  !$acc )
  do i = 1, n
    b(i)     = b(i)     * fn
    d_inv(i) = d_inv(i) * fn_inv
  end do
  !$acc end parallel

  end subroutine

