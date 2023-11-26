!==============================================================================!
  subroutine Write_Log(Comm, fh, var, disp)
!------------------------------------------------------------------------------!
!>  Writes a single logical variable to file in parallel environment.
!>  This subroutine is used only to write to backup files in Process and is
!>  invoked through Grid's member Comm in Backup_Mod.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type), intent(in)    :: Comm  !! communicator from grid
  type(Mpi_File),   intent(in)    :: fh    !! file handle
  logical,          intent(in)    :: var   !! variable to write
  integer(DP),      intent(inout) :: disp  !! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: error = 0
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Comm)
!==============================================================================!

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
                      var,                &
                      1,                  &
                      comm_type_log,      &
                      MPI_STATUS_IGNORE,  &
                      error)

  disp = disp + LP

  end subroutine
