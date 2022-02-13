!==============================================================================!
  subroutine Read_Real(Comm, fh, num, disp)
!------------------------------------------------------------------------------!
!   Read single real number for parallel runs.                                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type) :: Comm
  type(Mpi_File)   :: fh    ! file handle
  real             :: num   ! number to write out
  integer(DP)      :: disp  ! diplacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: error = 0
!==============================================================================!

  ! Set it at position disp (same as in Write counterpart)
  call Mpi_File_Set_View(fh,             &
                         disp,           &
                         MPI_DOUBLE,     &
                         MPI_DOUBLE,     &
                         'native',       &
                         MPI_INFO_NULL,  &
                         error)

  ! Read integer value 
  call Mpi_File_Read(fh,                 &
                     num,                &
                     1,                  &
                     MPI_DOUBLE,         &
                     MPI_STATUS_IGNORE,  &
                     error)

  disp = disp + SIZE_REAL

  end subroutine
