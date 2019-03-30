!==============================================================================!
  subroutine Comm_Mod_Write_Int_Array(fh, arr, disp)
!------------------------------------------------------------------------------!
!   Write integer array for parallel runs.                                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer               :: fh    ! file handle
  integer, dimension(:) :: arr   ! array to write out
  integer               :: disp  ! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: length
  integer :: error
!==============================================================================!

  ! Get array's length
  length = size(arr)

  ! Set it at position disp (same as in Read counterpart)
  call Mpi_File_Set_View(fh,             &
                         disp,           &
                         MPI_INTEGER8,   &
                         MPI_INTEGER8,   &
                         'native',       &
                         MPI_INFO_NULL,  &
                         error)

  ! Write integer value 
  call Mpi_File_Write(fh,                 &
                      arr(1),             &
                      length,             &
                      MPI_INTEGER8,       &
                      MPI_STATUS_IGNORE,  &
                      error) 

  disp = disp + SIZE_INT * length

  end subroutine
