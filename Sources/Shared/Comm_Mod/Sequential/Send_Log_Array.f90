!==============================================================================!
  subroutine Comm_Mod_Send_Log_Array(phi_s, len_s, dest)
!------------------------------------------------------------------------------!
!   Dummy function for sequential compilation.                                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical :: phi_s(len_s)  ! send buffer
  integer :: len_s         ! send length
  integer :: dest          ! destination processor
!==============================================================================!

  end subroutine
