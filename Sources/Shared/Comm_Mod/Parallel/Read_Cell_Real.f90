!==============================================================================!
  subroutine Read_Cell_Real(Comm, fh, array, disp)
!------------------------------------------------------------------------------!
!   Read distributed cell-based array.                                         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type) :: Comm
  type(Mpi_File)   :: fh         ! file handle
  real             :: array(:)
  integer(DP)      :: disp       ! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: error = 0
!==============================================================================!

  ! Set view for distributed cell data 
  ! (this part is the same as in Write counterpart)
  call Mpi_File_Set_View(fh,                    &
                         disp,                  &
                         MPI_DOUBLE,            &
                         Comm % cell_map_type,  &
                         'native',              &
                         MPI_INFO_NULL,         &
                         error)

  ! Read distributed cell data 
  call Mpi_File_Read(fh,                 &
                     array,              &
                     Comm % nc_sub,      &
                     MPI_DOUBLE,         &
                     MPI_STATUS_IGNORE,  &
                     error)

  disp = disp + Comm % nc_tot * SIZE_REAL

  end subroutine
