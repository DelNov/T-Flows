!==============================================================================!
  subroutine Exchange_Int_Array(Comm, length, phi, dest)
!------------------------------------------------------------------------------!
!   Dummy function for sequential runs.                                        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type) :: Comm
  integer          :: length
  integer          :: phi(length)
  integer          :: dest         ! destination processor
!==============================================================================!

  end subroutine
