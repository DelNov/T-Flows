!==============================================================================!
  subroutine Comm_Mod_Read_Real(fh, num, disp)
!------------------------------------------------------------------------------!
!   Read single real number for parallel runs.                                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: fh    ! file handle
  real    :: num   ! number to write out
  integer :: disp  ! diplacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: error
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
