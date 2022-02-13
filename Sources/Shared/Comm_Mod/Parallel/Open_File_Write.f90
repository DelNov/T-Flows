!==============================================================================!
  subroutine Open_File_Write(Comm, fh, file_name)
!------------------------------------------------------------------------------!
!   Open file for writing for parallel runs.                                   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type) :: Comm
  type(Mpi_File)   :: fh             ! file handle
  character        :: file_name*(*)  ! file name
!-----------------------------------[Locals]-----------------------------------!
  integer :: error = 0
!==============================================================================!

  ! Open file with MPI
  call Mpi_File_Open(MPI_COMM_WORLD,                    &
                     file_name,                         &
                     MPI_MODE_WRONLY + MPI_MODE_CREATE, & 
                     MPI_INFO_NULL,                     &
                     fh,                                &
                     error) 

  end subroutine
