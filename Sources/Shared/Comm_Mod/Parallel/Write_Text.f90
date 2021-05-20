!==============================================================================!
  subroutine Write_Text(Comm, fh, text_out, disp)
!------------------------------------------------------------------------------!
!   Write string array for parallel runs.                                      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type) :: Comm
  integer          :: fh            ! file handle
  character        :: text_out*(*)  ! text to write out
  integer          :: disp          ! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: leng, error
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
