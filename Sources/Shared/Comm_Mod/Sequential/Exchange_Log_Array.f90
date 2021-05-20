!==============================================================================!
  subroutine Exchange_Log_Array(Comm, length, phi, dest)
!------------------------------------------------------------------------------!
!   Dummy function for sequential runs.                                        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type) :: Comm
  integer          :: length
  logical          :: phi(length)
  integer          :: dest         ! destination processor
!==============================================================================!

  end subroutine
