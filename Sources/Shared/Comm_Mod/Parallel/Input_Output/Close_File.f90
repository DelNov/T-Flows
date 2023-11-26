!==============================================================================!
  subroutine Close_File(Comm, fh)
!------------------------------------------------------------------------------!
!>  Closes a file for parallel I/O. This subroutine is called from Process
!>  to close backup.  It is accessed through Grid's member Comm.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type), intent(in)    :: Comm  !! communicator from grid
  type(Mpi_File),   intent(inout) :: fh    !! file handle
!-----------------------------------[Locals]-----------------------------------!
  integer :: error
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Comm)
!==============================================================================!

  ! Close the file
  call Mpi_File_Close(fh, error)

  end subroutine
