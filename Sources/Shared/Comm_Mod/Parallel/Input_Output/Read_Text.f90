!==============================================================================!
  subroutine Read_Text(Comm, fh, text_in, disp)
!------------------------------------------------------------------------------!
!>  Reads a character array from a file in parallel environment. Each processor
!>  reads the same string, meaning that it will be distributed to all
!>  processors after reading. This subroutine is used only in Backup_Mod.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type), intent(in)    :: Comm         !! communicator from grid
  type(Mpi_File),   intent(in)    :: fh           !! file handle
  character,        intent(out)   :: text_in*(*)  !! text to read in
  integer(DP),      intent(inout) :: disp         !! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: leng
  integer :: error = 0
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Comm)
!==============================================================================!

  leng = len(text_in)

  ! Set it at position disp (same as in Write counterpart)
  call Mpi_File_Set_View(fh,             &
                         disp,           &
                         MPI_CHARACTER,  &
                         MPI_CHARACTER,  &
                         'native',       &
                         MPI_INFO_NULL,  &
                         error)

  ! Read the string
  call Mpi_File_Read(fh,                 &
                     text_in,            &
                     leng,               &
                     MPI_CHARACTER,      &
                     MPI_STATUS_IGNORE,  &
                     error) 

  disp = disp + leng

  end subroutine
