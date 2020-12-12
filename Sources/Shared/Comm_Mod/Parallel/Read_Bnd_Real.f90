!==============================================================================!
  subroutine Comm_Mod_Read_Bnd_Real(comm, fh, array, disp)
!------------------------------------------------------------------------------!
!   Read distributed boundary-cell-based array.                                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Comm_Type) :: comm
  integer         :: fh         ! file handle
  real            :: array(:)
  integer         :: disp       ! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: error
!==============================================================================!

  ! Set view for distributed boundary cell data 
  ! (this part is the same as in Write counterpart)
  call Mpi_File_Set_View(fh,                 &
                         disp,               &
                         MPI_DOUBLE,         &
                         bnd_cell_map_type,  &
                         'native',           &
                         MPI_INFO_NULL,      &
                         error)

  ! Read distributed boundary cell data 
  call Mpi_File_Read(fh,                 &
                     array,              &
                     comm % nb_sub,      &
                     MPI_DOUBLE,         &
                     MPI_STATUS_IGNORE,  &
                     error)

  disp = disp + comm % nb_tot * SIZE_REAL

  end subroutine
