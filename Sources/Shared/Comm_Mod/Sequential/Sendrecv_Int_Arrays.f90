!==============================================================================!
  subroutine Sendrecv_Int_Arrays(Comm, len_s, phi_s, &
                                       len_r, phi_r, dest)
!------------------------------------------------------------------------------!
!   Dummy function for sequential compilation.                                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type) :: Comm
  integer          :: len_s         ! send length
  integer          :: phi_s(len_s)  ! send buffer
  integer          :: len_r         ! receive length
  integer          :: phi_r(len_r)  ! receive buffer
  integer          :: dest          ! destination processor
!==============================================================================!

  end subroutine
