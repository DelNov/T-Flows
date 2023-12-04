!==============================================================================!
  logical function Probe_Int_Array(Mem, a, i)
!------------------------------------------------------------------------------!
!>  Checks if index i is within bounds of integer array a
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Memory_Type),   intent(in) :: Mem   !! parent class
  integer, allocatable, intent(in) :: a(:)  !! opearand array
  integer,              intent(in) :: i     !! array index
!==============================================================================!

  Probe_Int_Array = (i >= lbound(a, 1) .and. i <= ubound(a, 1))

  end function
