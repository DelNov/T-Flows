!==============================================================================!
  subroutine Send_Real_Array(Comm, len_s, phi_s, dest)
!------------------------------------------------------------------------------!
!   Dummy function for sequential compilation.                                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type) :: Comm
  integer          :: len_s         ! send length
  real             :: phi_s(len_s)  ! send buffer
  integer          :: dest          ! destination processor
!==============================================================================!

  end subroutine
