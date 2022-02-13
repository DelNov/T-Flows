!==============================================================================!
  subroutine Read_Real_Array(Comm, fh, arr, disp)
!------------------------------------------------------------------------------!
!   Read real array for parallel runs.                                         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type)     :: Comm
  type(Mpi_File)       :: fh    ! file handle
  real,   dimension(:) :: arr   ! array to read
  integer(DP)          :: disp  ! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: length
  integer :: error = 0
!==============================================================================!

  ! Get array's length
  length = size(arr)

  ! Set it at position disp (same as in Write counterpart)
  call Mpi_File_Set_View(fh,             &
                         disp,           &
                         MPI_DOUBLE,     &
                         MPI_DOUBLE,     &
                         'native',       &
                         MPI_INFO_NULL,  &
                         error)

  ! Read integer value 
  call Mpi_File_Read(fh,                 &
                     arr(1),             &
                     length,             &
                     MPI_DOUBLE,         &
                     MPI_STATUS_IGNORE,  &
                     error) 

  disp = disp + SIZE_REAL * length

  end subroutine
