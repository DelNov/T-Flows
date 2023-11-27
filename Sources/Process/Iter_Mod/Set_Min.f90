!==============================================================================!
  pure subroutine Set_Min(Iter, val)
!------------------------------------------------------------------------------!
!>  Function to set the minimum number of iterations in Iter.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Iter_Type), intent(inout) :: Iter  !! parent, singleton object Iter
  integer,          intent(in)    :: val   !! minimum number of iterations
!==============================================================================!

  Iter % min_iterations = val

  end subroutine
