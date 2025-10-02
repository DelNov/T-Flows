!==============================================================================!
  subroutine Write_Bnd_Real(Comm, fh, arr, disp)
!------------------------------------------------------------------------------!
!>  Writes a boundary-cell-based (hence, associated with a grid) real array
!>  from a parallel environment to a file. Each processor writes to a specific
!>  part of the file without interfering with others, which is achieved through
!>  field bnd_cell_map declared in Comm_Mod, but defined in the subroutine
!>  Grid % Form_Maps_For_Backup. This subroutine is used only to write to
!>  backup files and is invoked through the Grid's member Comm in Backup_Mod.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type), intent(in)    :: Comm    !! communicator from grid
  type(Mpi_File),   intent(in)    :: fh      !! file handle
  real,             intent(in)    :: arr(:)  !! bnd-cell-based array to write
  integer(DP),      intent(inout) :: disp    !! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: error = 0
!==============================================================================!

  ! Set view for distributed boundary cell data
  ! (this part is the same as in the Read counterpart)
  call Mpi_File_Set_View(fh,                        &
                         disp,                      &
                         comm_type_real,            &
                         Comm % bnd_cell_map_type,  &
                         'native',                  &
                         MPI_INFO_NULL,             &
                         error)

  ! Write distributed boundary cell data
  call Mpi_File_Write(fh,                 &
                      arr,                &
                      Comm % nb_sub,      &
                      comm_type_real,     &
                      MPI_STATUS_IGNORE,  &
                      error)

  disp = disp + Comm % nb_tot * RP

  end subroutine
