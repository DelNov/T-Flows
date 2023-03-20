!==============================================================================!
  subroutine Close_File(Comm, fh)
!------------------------------------------------------------------------------!
!   Close file for sequential runs.                                            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type) :: Comm
  integer          :: fh  ! file handle
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Comm)
!==============================================================================!

  close(fh)

  end subroutine
