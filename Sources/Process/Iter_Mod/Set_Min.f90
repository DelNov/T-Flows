!==============================================================================!
  pure subroutine Set_Min(Iter, val)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Iter_Type), intent(inout) :: Iter
  integer,          intent(in)    :: val
!==============================================================================!

  Iter % min_iterations = val

  end subroutine
