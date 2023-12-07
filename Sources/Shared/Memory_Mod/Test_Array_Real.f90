!==============================================================================!
  logical function Test_Array_Real(Mem, a, i)
!------------------------------------------------------------------------------!
!>  Checks if index i is within bounds of real array a
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Memory_Type), intent(in) :: Mem   !! parent class
  real, allocatable,  intent(in) :: a(:)  !! opearand array
  integer,            intent(in) :: i     !! array index
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Mem)
!==============================================================================!

  Test_Array_Real = (i >= lbound(a, 1) .and. i <= ubound(a, 1))

  end function
