!==============================================================================!
  subroutine Read_Log_Array(Comm, fh, arr, disp)
!------------------------------------------------------------------------------!
!   Read single integer for parallel runs.                                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type)      :: Comm
  integer               :: fh    ! file handle
  logical, dimension(:) :: arr   ! array to read
  integer               :: disp  ! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: length
  integer :: error
!==============================================================================!

  ! Get array's length
  length = size(arr)

  ! Set it at position disp (same as in Write counterpart)
  call Mpi_File_Set_View(fh,             &
                         disp,           &
                         MPI_LOGICAL,    &
                         MPI_LOGICAL,    &
                         'native',       &
                         MPI_INFO_NULL,  &
                         error)

  ! Read integer value 
  call Mpi_File_Read(fh,                 &
                     arr(1),             &
                     length,             &
                     MPI_LOGICAL,        &
                     MPI_STATUS_IGNORE,  &
                     error)

  disp = disp + SIZE_LOG * length

  end subroutine
