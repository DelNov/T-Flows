!==============================================================================!
  pure logical function First_Proc()
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  First_Proc = Global % this_processor < 2

  end function
