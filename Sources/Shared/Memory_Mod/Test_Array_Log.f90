!==============================================================================!
  logical function Test_Array_Log(Mem, a, i)
!------------------------------------------------------------------------------!
!>  Checks if index i is within bounds of logical array a
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Memory_Type),   intent(in) :: Mem   !! parent class
  logical, allocatable, intent(in) :: a(:)  !! opearand array
  integer,              intent(in) :: i     !! array index
!==============================================================================!

  Test_Array_Log = (i >= lbound(a, 1) .and. i <= ubound(a, 1))

  end function
