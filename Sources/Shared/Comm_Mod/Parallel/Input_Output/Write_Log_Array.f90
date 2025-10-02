!==============================================================================!
  subroutine Write_Log_Array(Comm, fh, arr, disp)
!------------------------------------------------------------------------------!
!>  Writes a logical array to file in parallel environment. Each processor
!>  writes the same data, so it better be distributed before writing.
!>  This subroutine is used only to write Swarm data to backup files and is
!>  invoked through Grid's member Comm in Backup_Mod.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type), intent(in)    :: Comm    !! communicator from grid
  type(Mpi_File),   intent(in)    :: fh      !! file handle
  logical,          intent(in)    :: arr(:)  !! array to write
  integer(DP),      intent(inout) :: disp    !! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: length
  integer :: error = 0
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Comm)
!==============================================================================!

  ! Get array's length
  length = size(arr)

  ! Set it at position disp (same as in Read counterpart)
  call Mpi_File_Set_View(fh,             &
                         disp,           &
                         comm_type_log,  &
                         comm_type_log,  &
                         'native',       &
                         MPI_INFO_NULL,  &
                         error)

  ! Write integer value 
  call Mpi_File_Write(fh,                 &
                      arr(1),             &
                      length,             &
                      comm_type_log,      &
                      MPI_STATUS_IGNORE,  &
                      error)

  disp = disp + LP * length

  end subroutine
