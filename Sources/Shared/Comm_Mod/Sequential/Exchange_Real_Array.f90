!==============================================================================!
  subroutine Exchange_Real_Array(Comm, length, phi, dest)
!------------------------------------------------------------------------------!
!   Dummy function for sequential compilation.                                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type) :: Comm
  integer          :: length
  real             :: phi(length)
  integer          :: dest         ! destination processor
!==============================================================================!

  end subroutine
