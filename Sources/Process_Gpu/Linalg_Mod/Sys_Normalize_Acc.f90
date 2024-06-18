!==============================================================================!
  subroutine Sys_Normalize_Acc(Lin, n, nz, fn, a_val, a_dia, d_inv, b)
!------------------------------------------------------------------------------!
!>  Calculation scalar over matrix diagonal operation on a device.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Linalg_Type) :: Lin         !! parent class
  integer            :: n           !! matrix and vector dimension
  integer            :: nz          !! number of nonzeros
  real,  intent(out) :: fn          !! resulting scaling factor
  real               :: a_val(nz)   !! operand matrix values
  integer            :: a_dia(n)    !! operand matrix positions of diagonals
  real               :: d_inv(n)    !! inverted diagonal entries
  real               :: b(n)        !! right hand side vector
!-----------------------------------[Locals]-----------------------------------!
  integer :: sum_n, i
  real    :: avg_a, sum_a, fn_inv
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Lin)
!==============================================================================!

  !-------------------------------------!
  !   Sum all the diagonal entries up   !
  !-------------------------------------!

  sum_a = 0.0

  !$acc parallel loop reduction(+:sum_a)
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

  !$acc parallel loop independent
  do i = 1, nz
    a_val(i) = a_val(i) * fn
  end do
  !$acc end parallel

  !$acc parallel loop independent
  do i = 1, n
    b(i)     = b(i)     * fn
    d_inv(i) = d_inv(i) * fn_inv
  end do
  !$acc end parallel

  end subroutine

