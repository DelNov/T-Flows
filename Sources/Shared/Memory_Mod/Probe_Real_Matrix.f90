!==============================================================================!
  logical function Probe_Real_Matrix(Mem, a, i, j)
!------------------------------------------------------------------------------!
!>  Checks if indices i and j are within bounds of real matrix a
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Memory_Type), intent(in) :: Mem     !! parent class
  real, allocatable,  intent(in) :: a(:,:)  !! operand matrix
  integer,            intent(in) :: i, j    !! matrix index
!==============================================================================!

  Probe_Real_Matrix = (i >= lbound(a, 1) .and. i <= ubound(a, 1) .and. &
                       j >= lbound(a, 2) .and. j <= ubound(a, 2))

  end function
