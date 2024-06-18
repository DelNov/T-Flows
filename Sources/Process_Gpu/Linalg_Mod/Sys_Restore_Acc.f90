!==============================================================================!
  subroutine Sys_Restore_Acc(Lin, n, nz, fn, a_val, a_dia, d_inv, b)
!------------------------------------------------------------------------------!
!>  Restores a linear system to its original (de-normalized) state.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Linalg_Type) :: Lin         !! parent class
  integer            :: n           !! matrix and vector dimension
  integer            :: nz          !! number of nonzeros
  real,  intent(in)  :: fn          !! resulting scaling factor
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

  !------------------------------------------------!
  !   Work out the inverse of the scaling factor   !
  !------------------------------------------------!
  fn_inv = 1.0 / fn       ! inverse scaling factor, like 1.0 / fn

  !-----------------------------------------------------!
  !   Scale the matrix and the right hand side vector   !
  !-----------------------------------------------------!

  !$acc parallel loop independent
  do i = 1, nz
    a_val(i) = a_val(i) * fn_inv
  end do
  !$acc end parallel

  !$acc parallel loop independent
  do i = 1, n
    b(i)     = b(i)     * fn_inv
    d_inv(i) = d_inv(i) * fn
  end do
  !$acc end parallel

  end subroutine

