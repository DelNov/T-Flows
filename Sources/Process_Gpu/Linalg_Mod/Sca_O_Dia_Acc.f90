!==============================================================================!
  subroutine Sca_O_Dia_Acc(Lin, n, nz, b, s, a_val, a_dia)
!------------------------------------------------------------------------------!
!>  Calculation scalar over matrix diagonal operation on a device.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Linalg_Type) :: Lin         !! parent class
  integer            :: n           !! matrix and vector dimension
  integer            :: nz          !! number of nonzeros
  real               :: b(n)        !! result vector
  real               :: s           !! scalar
  real               :: a_val(nz)   !! operand matrix values
  integer            :: a_dia(n)    !! operand matrix positions of diagonals
!-----------------------------------[Locals]-----------------------------------!
  integer :: i
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Lin)
!==============================================================================!

  !$acc kernels present(b, a_val, a_dia)
  do i = 1, n
    b(i) = s / a_val(a_dia(i))
  end do
  !$acc end kernels

  end subroutine

