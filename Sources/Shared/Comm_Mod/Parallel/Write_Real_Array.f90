!==============================================================================!
  subroutine Write_Real_Array(Comm, fh, arr, disp)
!------------------------------------------------------------------------------!
!   Write real arrayr for parallel runs.                                       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type),   intent(in)    :: Comm
  type(Mpi_File),     intent(in)    :: fh    ! file handle
  real, dimension(:), intent(in)    :: arr   ! array to write out
  integer(DP),        intent(inout) :: disp  ! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: length
  integer :: error = 0
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Comm)
!==============================================================================!

  ! Get array's length
  length = size(arr)

  ! Set it at position disp (same as in Read counterpart)
  call Mpi_File_Set_View(fh,             &
                         disp,           &
                         comm_type_real, &
                         comm_type_real, &
                         'native',       &
                         MPI_INFO_NULL,  &
                         error)

  ! Write real value 
  call Mpi_File_Write(fh,                 &
                      arr(1),             &
                      length,             &
                      comm_type_real,     &
                      MPI_STATUS_IGNORE,  &
                      error)

  disp = disp + RP * length

  end subroutine
