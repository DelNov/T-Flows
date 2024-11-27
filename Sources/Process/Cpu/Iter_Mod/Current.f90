!==============================================================================!
  pure integer function Current(Iter)
!------------------------------------------------------------------------------!
!>  Function to retreive the current iteration stored in Iter.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Iter_Type), intent(in) :: Iter  !! parent, singleton object Iter
!==============================================================================!

  Current = Iter % current_iteration

  end function
