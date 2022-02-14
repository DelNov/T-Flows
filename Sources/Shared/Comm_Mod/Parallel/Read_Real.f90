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
                         comm_type_real, &
                         comm_type_real, &
                         'native',       &
                         MPI_INFO_NULL,  &
                         error)

  ! Read integer value 
  call Mpi_File_Read(fh,                 &
                     num,                &
                     1,                  &
                     comm_type_real,     &
                     MPI_STATUS_IGNORE,  &
                     error)

  disp = disp + RP

  end subroutine
