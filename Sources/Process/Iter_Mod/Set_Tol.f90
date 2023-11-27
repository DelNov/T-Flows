!==============================================================================!
  pure subroutine Set_Tol(Iter, val)
!------------------------------------------------------------------------------!
!>  Function to set the iteration tolerance (for SIMPLE/PISO algorithm).
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Iter_Type), intent(inout) :: Iter  !! parent, singleton object Iter
  real,             intent(in)    :: val   !! desired tolerance
!==============================================================================!

  Iter % tolerance = val

  end subroutine
