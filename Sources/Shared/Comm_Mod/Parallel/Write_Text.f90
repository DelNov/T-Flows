!==============================================================================!
  subroutine Write_Text(Comm, fh, text_out, disp)
!------------------------------------------------------------------------------!
!   Write string array for parallel runs.                                      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type) :: Comm
  type(Mpi_File)   :: fh            ! file handle
  character        :: text_out*(*)  ! text to write out
  integer(DP)      :: disp          ! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: leng
  integer :: error = 0
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Comm)
!==============================================================================!

  leng = len(text_out)

  ! Set it at position disp (same as in Read counterpart)
  call Mpi_File_Set_View(fh,             &
                         disp,           &
                         MPI_CHARACTER,  &
                         MPI_CHARACTER,  &
                         'native',       &
                         MPI_INFO_NULL,  &
                         error)

  ! Write string out
  call Mpi_File_Write(fh,                 &
                      text_out,           &
                      leng,               &
                      MPI_CHARACTER,      &
                      MPI_STATUS_IGNORE,  &
                      error)

  disp = disp + leng

  end subroutine
