!==============================================================================!
  subroutine Read_Log(Comm, fh, var, disp)
!------------------------------------------------------------------------------!
!   Read single logical variable for parallel runs                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type) :: Comm
  type(Mpi_File)   :: fh    ! file handle
  logical          :: var   ! variable to read
  integer(DP)      :: disp  ! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: error = 0
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
