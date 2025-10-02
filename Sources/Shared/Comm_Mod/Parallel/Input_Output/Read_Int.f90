!==============================================================================!
  subroutine Read_Int(Comm, fh, num, disp)
!------------------------------------------------------------------------------!
!>  Reads an integer from a file in parallel environment. Each processor
!>  reads the same number, meaning that it will be distributed to all
!>  processors after reading. This subroutine is used only in Backup_Mod.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type), intent(in)    :: Comm  !! communicator from grid
  type(Mpi_File),   intent(in)    :: fh    !! file handle
  integer,          intent(out)   :: num   !! number to read
  integer(DP),      intent(inout) :: disp  !! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: error = 0
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Comm)
!==============================================================================!

  ! Set it at position disp (same as in Write counterpart)
  call Mpi_File_Set_View(fh,             &
                         disp,           &
                         comm_type_int,  &
                         comm_type_int,  &
                         'native',       &
                         MPI_INFO_NULL,  &
                         error)

  ! Read integer value
  call Mpi_File_Read(fh,                 &
                     num,                &
                     1,                  &
                     comm_type_int,      &
                     MPI_STATUS_IGNORE,  &
                     error)

  disp = disp + IP

  end subroutine
