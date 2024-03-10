!==============================================================================!
  pure subroutine Set_Max(Iter, val)
!------------------------------------------------------------------------------!
!>  Function to set the maximum number of iterations in Iter.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Iter_Type), intent(inout) :: Iter  !! parent, singleton object Iter
  integer,          intent(in)    :: val   !! maximum number of iterations
!==============================================================================!

  Iter % max_iterations = val

  end subroutine
