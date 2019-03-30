!==============================================================================!
  subroutine Comm_Mod_Write_Log_Array(fh, arr, disp)
!------------------------------------------------------------------------------!
!   Write integer array for parallel runs.                                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer               :: fh    ! file handle
  logical, dimension(:) :: arr   ! array to write out
  integer               :: disp  ! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: length
  integer :: error
!==============================================================================!

  ! Get array's length
  length = size(arr)

  ! Set it at position disp (same as in Read counterpart)
  call Mpi_File_Set_View(fh,             &
                         disp,           &
                         MPI_LOGICAL8,   &
                         MPI_LOGICAL8,   &
                         'native',       &
                         MPI_INFO_NULL,  &
                         error)

  ! Write integer value 
  call Mpi_File_Write(fh,                 &
                      arr(1),             &
                      length,             &
                      MPI_LOGICAL8,       &
                      MPI_STATUS_IGNORE,  &
                      error) 

  disp = disp + SIZE_LOG * length

  end subroutine
