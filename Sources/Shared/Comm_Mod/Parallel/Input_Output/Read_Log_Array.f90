!==============================================================================!
  subroutine Read_Log_Array(Comm, fh, arr, disp)
!------------------------------------------------------------------------------!
!>  Reads a logical array from a file in parallel environment. Each processor
!>  reads the entire array, meaning that the data will be distributed to all
!>  processors after reading. This subroutine is used only to read Swarm data
!>  from backup files and is invoked through Grid's member Comm in Backup_Mod.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type), intent(in)    :: Comm    !! communicator from grid
  type(Mpi_File),   intent(in)    :: fh      !! file handle
  logical,          intent(out)   :: arr(:)  !! array to read
  integer(DP),      intent(inout) :: disp    !! displacement in bytes
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
                         comm_type_log,  &
                         comm_type_log,  &
                         'native',       &
                         MPI_INFO_NULL,  &
                         error)

  ! Read integer value
  call Mpi_File_Read(fh,                 &
                     arr(1),             &
                     length,             &
                     comm_type_log,      &
                     MPI_STATUS_IGNORE,  &
                     error)

  disp = disp + LP * length

  end subroutine
