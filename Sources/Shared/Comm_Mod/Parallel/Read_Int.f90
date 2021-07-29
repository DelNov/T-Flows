!==============================================================================!
  subroutine Read_Int(Comm, fh, num, disp)
!------------------------------------------------------------------------------!
!   Read single integer for parallel runs.                                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type) :: Comm
  integer          :: fh    ! file handle
  integer          :: num   ! number to read
  integer          :: disp  ! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: error = 0
!==============================================================================!

  ! Set it at position disp (same as in Write counterpart)
  call Mpi_File_Set_View(fh,             &   
                         disp,           &   
                         MPI_INTEGER8,   &
                         MPI_INTEGER8,   &   
                         'native',       &   
                         MPI_INFO_NULL,  &
                         error)

  ! Read integer value 
  call Mpi_File_Read(fh,                 &
                     num,                &
                     1,                  &
                     MPI_INTEGER8,       &
                     MPI_STATUS_IGNORE,  &
                     error) 

  disp = disp + SIZE_INT

  end subroutine
