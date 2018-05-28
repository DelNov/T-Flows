!==============================================================================!
  subroutine Comm_Mod_Write_Bnd_Real(fh, array, disp)
!------------------------------------------------------------------------------!
!   Write distributed boundary-cell-based array.                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: fh         ! file handle
  real    :: array(:)
  integer :: disp       ! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: error
!==============================================================================!

  ! Set view for distributed boundary cell data 
  ! (this part is the same as in the Read counterpart)
  call Mpi_File_Set_View(fh,                 &   
                         disp,               &   
                         MPI_DOUBLE,         &
                         bnd_cell_map_type,  &   
                         'native',           &   
                         MPI_INFO_NULL,      &
                         error)

  ! Write distributed boundary cell data 
  call Mpi_File_Write(fh,                 &
                      array,              &
                      nb_s,               &   
                      MPI_DOUBLE,         &   
                      MPI_STATUS_IGNORE,  &
                      error)

  disp = disp + nb_t * SIZE_REAL

  end subroutine
