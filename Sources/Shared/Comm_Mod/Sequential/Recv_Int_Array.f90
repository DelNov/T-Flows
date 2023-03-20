!==============================================================================!
  subroutine Recv_Int_Array(Comm, len_r, phi_r, dest)
!------------------------------------------------------------------------------!
!   Dummy function for sequential compilation.                                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type) :: Comm
  integer          :: len_r         ! receive length
  integer          :: phi_r(len_r)  ! receive buffer
  integer          :: dest          ! destination processor
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Comm)
  Unused(len_r)
  Unused(phi_r)
  Unused(dest)
!==============================================================================!

  end subroutine
