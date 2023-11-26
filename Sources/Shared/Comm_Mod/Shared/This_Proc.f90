!==============================================================================!
  pure integer function This_Proc()
!------------------------------------------------------------------------------!
!>  This utility function is used for identifying the rank of the current
!>  processor in a parallel computing environment. It achieves it by returning
!>  the value of the private member this_processor from the global Comm_Type
!>  object, Global. This function is crucial in scenarios where specific
!>  actions or decisions need to be made based on the rank of the processor
!>  executing the code, ensuring proper coordination among processors in a
!>  parallel computation setting.
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  This_Proc = Global % this_processor

  end function
