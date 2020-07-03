!==============================================================================!
  subroutine Comm_Mod_Send_Real_Array(phi_s, len_s, dest)
!------------------------------------------------------------------------------!
!   Dummy function for sequential compilation.                                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real    :: phi_s(len_s)  ! send buffer
  integer :: len_s         ! send length
  integer :: dest          ! destination processor
!==============================================================================!

  end subroutine
