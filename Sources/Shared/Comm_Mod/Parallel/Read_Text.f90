!==============================================================================!
  subroutine Comm_Mod_Read_Text(fh, text_in, disp)
!------------------------------------------------------------------------------!
!   Read string array for paralle runs.                                        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer   :: fh           ! file handle
  character :: text_in*(*)  ! text to write out
  integer   :: disp         ! displacement
!-----------------------------------[Locals]-----------------------------------!
  integer :: leng, error
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
