!==============================================================================!
  pure integer function Get_Max(Iter)
!------------------------------------------------------------------------------!
!>  Function to retreive the maximum number of iterations from Iter.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Iter_Type), intent(in) :: Iter  !! parent, singleton object Iter
!==============================================================================!

  Get_Max = Iter % max_iterations

  end function
