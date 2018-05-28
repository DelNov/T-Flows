!==============================================================================!
  subroutine Comm_Mod_Open_File_Read(fh, file_name)
!------------------------------------------------------------------------------!
!   Open file for reading for parallel runs.                                   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer   :: fh             ! file handle
  character :: file_name*(*)  ! file name
!-----------------------------------[Locals]-----------------------------------!
  integer :: error
!==============================================================================!

  ! Open file with MPI
  call Mpi_File_Open(MPI_COMM_WORLD,   &
                     file_name,        &
                     MPI_MODE_RDONLY,  &
                     MPI_INFO_NULL,    &
                     fh,               &
                     error) 

  end subroutine
