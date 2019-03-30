!==============================================================================!
  subroutine Comm_Mod_Write_Real_Array(fh, arr, disp)
!------------------------------------------------------------------------------!
!   Write real arrayr for parallel runs.                                       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer              :: fh    ! file handle
  real,   dimension(:) :: arr   ! array to write out
  integer              :: disp  ! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: length
  integer :: error
!==============================================================================!

  ! Get array's length
  length = size(arr)

  ! Set it at position disp (same as in Read counterpart)
  call Mpi_File_Set_View(fh,             &
                         disp,           &
                         MPI_DOUBLE,     &
                         MPI_DOUBLE,     &
                         'native',       &
                         MPI_INFO_NULL,  &
                         error)

  ! Write real value 
  call Mpi_File_Write(fh,                 &
                      arr(1),             &
                      length,             &
                      MPI_DOUBLE,         &
                      MPI_STATUS_IGNORE,  &
                      error) 

  disp = disp + SIZE_REAL * length

  end subroutine
