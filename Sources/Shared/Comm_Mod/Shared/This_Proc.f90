!==============================================================================!
  pure integer function This_Proc()
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  This_Proc = Global % this_processor

  end function
