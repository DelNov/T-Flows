!==============================================================================!
  subroutine Open_File_Read(Comm, fh, file_name)
!------------------------------------------------------------------------------!
!>  Opens a file for reading for sequential I/O. This subroutine is called
!>  from Process to open backup files for reading.  It is accessed through
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

  open(newunit=fh, FILE=file_name, access='stream', form='unformatted')

  end subroutine
