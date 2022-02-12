!==============================================================================!
  subroutine Read_Text(Comm, fh, text_in, disp)
!------------------------------------------------------------------------------!
!   Read string array for paralle runs.                                        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type) :: Comm
  integer          :: fh           ! file handle
  character        :: text_in*(*)  ! text to write out
  integer(DP)      :: disp         ! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: leng
  integer :: error = 0
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
