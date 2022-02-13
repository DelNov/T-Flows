!==============================================================================!
  subroutine Write_Int(Comm, fh, num, disp)
!------------------------------------------------------------------------------!
!   Write single integer for parallel runs.                                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type) :: Comm
  type(Mpi_File)   :: fh    ! file handle
  integer          :: num   ! number to write out
  integer(DP)      :: disp  ! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: error = 0
!==============================================================================!

  ! Set it at position disp (same as in Read counterpart)
  call Mpi_File_Set_View(fh,             &
                         disp,           &
                         MPI_INTEGER,    &
                         MPI_INTEGER,    &
                         'native',       &
                         MPI_INFO_NULL,  &
                         error)

  ! Write integer value 
  call Mpi_File_Write(fh,                 &
                      num,                &
                      1,                  &
                      MPI_INTEGER,        &
                      MPI_STATUS_IGNORE,  &
                      error) 

  disp = disp + SIZE_INT

  end subroutine
