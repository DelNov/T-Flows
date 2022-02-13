!==============================================================================!
  subroutine Read_Int_Array(Comm, fh, arr, disp)
!------------------------------------------------------------------------------!
!   Read integer array for parallel runs                                       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type)      :: Comm
  type(Mpi_File)        :: fh    ! file handle
  integer, dimension(:) :: arr   ! array to read
  integer(DP)           :: disp  ! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: length
  integer :: error = 0
!==============================================================================!

  ! Get array's length
  length = size(arr)

  ! Set it at position disp (same as in Write counterpart)
  call Mpi_File_Set_View(fh,             &
                         disp,           &
                         MPI_INTEGER,    &
                         MPI_INTEGER,    &
                         'native',       &
                         MPI_INFO_NULL,  &
                         error)

  ! Read integer value 
  call Mpi_File_Read(fh,                 &
                     arr(1),             &
                     length,             &
                     MPI_INTEGER,        &
                     MPI_STATUS_IGNORE,  &
                     error)

  disp = disp + SIZE_INT * length

  end subroutine
