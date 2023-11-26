!==============================================================================!
  subroutine Open_File_Write(Comm, fh, file_name)
!------------------------------------------------------------------------------!
!>  Opens a file for writing for sequential I/O. This subroutine is called
!>  from Process to open backup files for writing.  It is accessed through
!>  Grid's member Comm.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type), intent(in)  :: Comm           !! communicator from grid
  integer,          intent(out) :: fh             !! file handle
  character,        intent(in)  :: file_name*(*)  !! file name
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Comm)
!==============================================================================!

  open(newunit=fh, file=file_name, access='stream', form='unformatted')

  end subroutine
