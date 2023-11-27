!==============================================================================!
  pure integer function Get_Min(Iter)
!------------------------------------------------------------------------------!
!>  Function to retreive the minimum number of iterations from Iter.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Iter_Type), intent(in) :: Iter  !! parent, singleton object Iter
!==============================================================================!

  Get_Min = Iter % min_iterations

  end function
