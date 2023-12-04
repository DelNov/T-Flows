!==============================================================================!
  logical function Probe_Real_Array(Mem, a, i)
!------------------------------------------------------------------------------!
!>  Checks if index i is within bounds of real array a
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Memory_Type), intent(in) :: Mem   !! parent class
  real, allocatable,  intent(in) :: a(:)  !! opearand array
  integer,            intent(in) :: i     !! array index
!==============================================================================!

  Probe_Real_Array = (i >= lbound(a, 1) .and. i <= ubound(a, 1))

  end function
