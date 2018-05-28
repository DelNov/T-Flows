!==============================================================================!
  subroutine Comm_Mod_Write_Real(fh, num, disp)
!------------------------------------------------------------------------------!
!   Write single reael number for parallel runs.                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: fh    ! file handle
  real    :: num   ! number to write out
  integer :: disp  ! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: error
!==============================================================================!

  ! Set it at position disp (same as in Read counterpart)
  call Mpi_File_Set_View(fh,             &   
                         disp,           &   
                         MPI_DOUBLE,     &
                         MPI_DOUBLE,     &   
                         'native',       &   
                         MPI_INFO_NULL,  &
                         error)

  ! Write integer value 
  call Mpi_File_Write(fh,                 &
                      num,                &
                      1,                  &
                      MPI_DOUBLE,         &
                      MPI_STATUS_IGNORE,  &
                      error) 

  disp = disp + SIZE_REAL

  end subroutine
