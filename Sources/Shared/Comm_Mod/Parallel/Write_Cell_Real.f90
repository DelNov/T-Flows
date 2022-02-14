!==============================================================================!
  subroutine Write_Cell_Real(Comm, fh, array, disp)
!------------------------------------------------------------------------------!
!   Write distributed cell-based array.                                        !
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
  ! (this part is the same as in Read counterpart)
  call Mpi_File_Set_View(fh,                    &
                         disp,                  &
                         comm_type_real,        &
                         Comm % cell_map_type,  &
                         'native',              &
                         MPI_INFO_NULL,         &
                         error)

  ! Write distributed cell data 
  call Mpi_File_Write(fh,                 &
                      array,              &
                      Comm % nc_sub,      &
                      comm_type_real,     &
                      MPI_STATUS_IGNORE,  &
                      error)

  disp = disp + Comm % nc_tot * RP

  end subroutine
