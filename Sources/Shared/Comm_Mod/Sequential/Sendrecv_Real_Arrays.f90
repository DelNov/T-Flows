!==============================================================================!
  subroutine Sendrecv_Real_Arrays(Comm, len_s, phi_s,  &
                                        len_r, phi_r, dest)
!------------------------------------------------------------------------------!
!   Dummy function for sequential compilation.                                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type) :: Comm
  integer          :: len_s         ! send length
  real             :: phi_s(len_s)  ! send buffer
  integer          :: len_r         ! receive length
  real             :: phi_r(len_r)  ! receive buffer
  integer          :: dest          ! destination processor
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Comm)
  Unused(len_s)
  Unused(phi_s)
  Unused(len_r)
  Unused(phi_r)
  Unused(dest)
!==============================================================================!

  end subroutine
