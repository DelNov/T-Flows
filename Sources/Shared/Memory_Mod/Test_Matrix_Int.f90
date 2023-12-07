!==============================================================================!
  logical function Test_Matrix_Int(Mem, a, i, j)
!------------------------------------------------------------------------------!
!>  Checks if indices i and j are within bounds of integer matrix a
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Memory_Type),   intent(in) :: Mem     !! parent class
  integer, allocatable, intent(in) :: a(:,:)  !! operand matrix
  integer,              intent(in) :: i, j    !! matrix index
!==============================================================================!

  Test_Matrix_Int = (i >= lbound(a, 1) .and. i <= ubound(a, 1) .and. &
                     j >= lbound(a, 2) .and. j <= ubound(a, 2))

  end function
