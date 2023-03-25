!==============================================================================!
  subroutine Read_Cell_Real(Comm, fh, array, disp)
!------------------------------------------------------------------------------!
!   Read distributed cell-based array.                                         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type),   intent(in)    :: Comm
  type(Mpi_File),     intent(in)    :: fh         ! file handle
  real, dimension(:), intent(out)   :: array(:)   ! array to read
  integer(DP),        intent(inout) :: disp       ! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: error = 0
!==============================================================================!

  ! Set view for distributed cell data 
  ! (this part is the same as in Write counterpart)
  call Mpi_File_Set_View(fh,                    &
                         disp,                  &
                         comm_type_real,        &
                         Comm % cell_map_type,  &
                         'native',              &
                         MPI_INFO_NULL,         &
                         error)

  ! Read distributed cell data 
  call Mpi_File_Read(fh,                 &
                     array,              &
                     Comm % nc_sub,      &
                     comm_type_real,     &
                     MPI_STATUS_IGNORE,  &
                     error)

  disp = disp + Comm % nc_tot * RP

  end subroutine
