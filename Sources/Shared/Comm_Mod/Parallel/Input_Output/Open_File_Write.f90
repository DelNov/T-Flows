!==============================================================================!
  subroutine Open_File_Write(Comm, fh, file_name)
!------------------------------------------------------------------------------!
!>  Opens a file for writing for parallel I/O. This subroutine is called from
!>  Process to open backup files for writing.  It is accessed through Grid's
!>  member Comm.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type), intent(in)  :: Comm           !! communicator from grid
  type(Mpi_File),   intent(out) :: fh             !! file handle
  character,        intent(in)  :: file_name*(*)  !! file name
!-----------------------------------[Locals]-----------------------------------!
  integer :: error = 0
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Comm)
!==============================================================================!

  ! Open file with MPI
  call Mpi_File_Open(MPI_COMM_WORLD,                    &
                     file_name,                         &
                     MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                     MPI_INFO_NULL,                     &
                     fh,                                &
                     error)

  end subroutine
