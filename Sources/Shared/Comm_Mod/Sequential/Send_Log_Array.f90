!==============================================================================!
  subroutine Send_Log_Array(Comm, len_s, phi_s, dest)
!------------------------------------------------------------------------------!
!   Dummy function for sequential compilation.                                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type) :: Comm
  integer          :: len_s         ! send length
  logical          :: phi_s(len_s)  ! send buffer
  integer          :: dest          ! destination processor
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Comm)
  Unused(len_s)
  Unused(phi_s)
  Unused(dest)
!==============================================================================!

  end subroutine
