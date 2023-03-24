!==============================================================================!
  subroutine Sendrecv_Log_Arrays(Comm, len_s, phi_s,  &
                                       len_r, phi_r, dest)
!------------------------------------------------------------------------------!
!   Dummy function for sequential compilation.                                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type), intent(in)  :: Comm
  integer,          intent(in)  :: len_s         ! send length
  logical,          intent(in)  :: phi_s(len_s)  ! send buffer
  integer,          intent(in)  :: len_r         ! receive length
  logical,          intent(out) :: phi_r(len_r)  ! receive buffer
  integer,          intent(in)  :: dest          ! destination processor
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Comm)
  Unused(len_s)
  Unused(phi_s)
  Unused(len_r)
  Unused(phi_r)
  Unused(dest)
!==============================================================================!

  end subroutine
