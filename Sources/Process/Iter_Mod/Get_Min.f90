!==============================================================================!
  pure integer function Get_Min(Iter)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Iter_Type), intent(in) :: Iter
!==============================================================================!

  Get_Min = Iter % min_iterations

  end function
