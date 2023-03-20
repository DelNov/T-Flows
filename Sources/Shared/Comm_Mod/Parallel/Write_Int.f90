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
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Comm)
!==============================================================================!

  ! Set it at position disp (same as in Read counterpart)
  call Mpi_File_Set_View(fh,             &
                         disp,           &
                         comm_type_int,  &
                         comm_type_int,  &
                         'native',       &
                         MPI_INFO_NULL,  &
                         error)

  ! Write integer value 
  call Mpi_File_Write(fh,                 &
                      num,                &
                      1,                  &
                      comm_type_int,      &
                      MPI_STATUS_IGNORE,  &
                      error) 

  disp = disp + IP

  end subroutine
