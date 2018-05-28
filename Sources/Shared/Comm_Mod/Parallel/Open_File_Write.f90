!==============================================================================!
  subroutine Comm_Mod_Open_File_Write(fh, file_name)
!------------------------------------------------------------------------------!
!   Open file for writing for parallel runs.                                   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer   :: fh             ! file handle
  character :: file_name*(*)  ! file name
!-----------------------------------[Locals]-----------------------------------!
  integer :: error
!==============================================================================!

  ! Open file with MPI
  call Mpi_File_Open(MPI_COMM_WORLD,                    &
                     file_name,                         &
                     MPI_MODE_WRONLY + MPI_MODE_CREATE, & 
                     MPI_INFO_NULL,                     &
                     fh,                                &
                     error) 

  end subroutine
