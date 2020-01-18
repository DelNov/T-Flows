!==============================================================================!
  subroutine Comm_Mod_Write_Cell_Real(comm, fh, array, disp)
!------------------------------------------------------------------------------!
!   Write distributed cell-based array.                                        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Comm_Type) :: comm
  integer         :: fh         ! file handle
  real            :: array(:)
  integer         :: disp       ! displacement
!-----------------------------------[Locals]-----------------------------------!
  integer :: error
!==============================================================================!

  ! Set view for distributed cell data 
  ! (this part is the same as in Read counterpart)
  call Mpi_File_Set_View(fh,             &
                         disp,           &
                         MPI_DOUBLE,     &
                         cell_map_type,  &
                         'native',       &
                         MPI_INFO_NULL,  &
                         error)

  ! Write distributed cell data 
  call Mpi_File_Write(fh,                 &
                      array,              &
                      comm % nc_s,        &
                      MPI_DOUBLE,         &
                      MPI_STATUS_IGNORE,  &
                      error)

  disp = disp + comm % nc_t * SIZE_REAL

  end subroutine
