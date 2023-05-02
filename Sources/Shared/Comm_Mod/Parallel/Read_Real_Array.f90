!==============================================================================!
  subroutine Read_Real_Array(Comm, fh, arr, disp)
!------------------------------------------------------------------------------!
!   Read real array for parallel runs.                                         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type),   intent(in)    :: Comm
  type(Mpi_File),     intent(in)    :: fh    ! file handle
  real, dimension(:), intent(out)   :: arr   ! array to read
  integer(DP),        intent(inout) :: disp  ! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: length
  integer :: error = 0
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Comm)
!==============================================================================!

  ! Get array's length
  length = size(arr)

  ! Set it at position disp (same as in Write counterpart)
  call Mpi_File_Set_View(fh,             &
                         disp,           &
                         comm_type_real, &
                         comm_type_real, &
                         'native',       &
                         MPI_INFO_NULL,  &
                         error)

  ! Read integer value 
  call Mpi_File_Read(fh,                 &
                     arr(1),             &
                     length,             &
                     comm_type_real,     &
                     MPI_STATUS_IGNORE,  &
                     error) 

  disp = disp + RP * length

  end subroutine
