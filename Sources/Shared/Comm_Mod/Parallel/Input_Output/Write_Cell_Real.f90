!==============================================================================!
  subroutine Write_Cell_Real(Comm, fh, arr, disp)
!------------------------------------------------------------------------------!
!>  Writes a cell-based (hence, associated with a grid) real array from a
!>  parallel environment to a file. Each processor writes to a specific part
!>  of the file without interfering with others, which is achieved through
!>  fields cell_map declared in Comm_Mod, but defined in the subroutine
!>  Grid % Form_Maps_For_Backup. This subroutine is used only to write to
!>  backup files and is invoked through the Grid's member Comm in Backup_Mod.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type), intent(in)    :: Comm    !! communicator from the grid
  type(Mpi_File),   intent(in)    :: fh      !! file handle
  real,             intent(out)   :: arr(:)  !! bnd-cell-based array to write
  integer(DP),      intent(inout) :: disp    !! displacement in bytes
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
                      arr,                &
                      Comm % nc_sub,      &
                      comm_type_real,     &
                      MPI_STATUS_IGNORE,  &
                      error)

  disp = disp + Comm % nc_tot * RP

  end subroutine
