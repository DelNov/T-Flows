!==============================================================================!
  subroutine Comm_Mod_Sendrecv_Int_Arrays(len_s, phi_s, &
                                          len_r, phi_r, dest)
!------------------------------------------------------------------------------!
!   Dummy function for sequential compilation.                                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: len_s         ! send length
  integer :: phi_s(len_s)  ! send buffer
  integer :: len_r         ! receive length
  integer :: phi_r(len_r)  ! receive buffer
  integer :: dest          ! destination processor
!==============================================================================!

  end subroutine
