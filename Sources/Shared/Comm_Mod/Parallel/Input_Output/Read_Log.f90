!==============================================================================!
  subroutine Read_Log(Comm, fh, var, disp)
!------------------------------------------------------------------------------!
!>  Reads a logical value from a file in parallel environment. Each processor
!>  reads the same value, meaning that it will be distributed to all
!>  processors after reading. This subroutine is used only in Backup_Mod.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type), intent(in)    :: Comm  !! communicator from grid
  type(Mpi_File),   intent(in)    :: fh    !! file handle
  logical,          intent(out)   :: var   !! variable to read
  integer(DP),      intent(inout) :: disp  !! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: error = 0
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Comm)
!==============================================================================!

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
                     var,                &
                     1,                  &
                     comm_type_log,      &
                     MPI_STATUS_IGNORE,  &
                     error)

  disp = disp + LP

  end subroutine
