!==============================================================================!
  subroutine Write_Int(Comm, fh, num, disp)
!------------------------------------------------------------------------------!
!>  Writes a single integer number to a file in parallel environment.
!>  This subroutine is used only to write to backup files in Process and is
!>  invoked through Grid's member Comm in Backup_Mod.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type), intent(in)    :: Comm  !! communicator from grid
  type(Mpi_File),   intent(in)    :: fh    !! file handle
  integer,          intent(in)    :: num   !! number to write
  integer(DP),      intent(inout) :: disp  !! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: error = 0
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Comm)
!==============================================================================!

  ! Set it at position disp (same as in Read counterpart)
  call Mpi_File_Set_View(fh,             &
                         disp,           &
                         comm_type_int,  &
                         comm_type_int,  &
                         'native',       &
                         MPI_INFO_NULL,  &
                         error)

  ! Write integer value 
  call Mpi_File_Write(fh,                 &
                      num,                &
                      1,                  &
                      comm_type_int,      &
                      MPI_STATUS_IGNORE,  &
                      error) 

  disp = disp + IP

  end subroutine
