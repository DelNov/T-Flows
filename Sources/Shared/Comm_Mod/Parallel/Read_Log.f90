!==============================================================================!
  subroutine Read_Log(Comm, fh, var, disp)
!------------------------------------------------------------------------------!
!   Read single logical variable for parallel runs                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type) :: Comm
  integer          :: fh    ! file handle
  logical          :: var   ! variable to read
  integer          :: disp  ! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: error = 0
!==============================================================================!

  ! Set it at position disp (same as in Write counterpart)
  call Mpi_File_Set_View(fh,             &
                         disp,           &
                         MPI_LOGICAL8,   &
                         MPI_LOGICAL8,   &
                         'native',       &
                         MPI_INFO_NULL,  &
                         error)

  ! Read integer value 
  call Mpi_File_Read(fh,                 &
                     var,                &
                     1,                  &
                     MPI_LOGICAL8,       &
                     MPI_STATUS_IGNORE,  &
                     error)

  disp = disp + SIZE_LOG

  end subroutine
