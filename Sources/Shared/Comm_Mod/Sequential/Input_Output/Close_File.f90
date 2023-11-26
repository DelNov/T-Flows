!==============================================================================!
  subroutine Close_File(Comm, fh)
!------------------------------------------------------------------------------!
!>  Closes a file for sequential I/O. This subroutine is called from Process
!>  to close backup.  It is accessed through Grid's member Comm.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type), intent(in)    :: Comm  !! communicator from grid
  integer,          intent(inout) :: fh    !! file handle
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Comm)
!==============================================================================!

  close(fh)

  end subroutine
