!==============================================================================!
  subroutine Write_Bnd_Real(Comm, fh, array, disp)
!------------------------------------------------------------------------------!
!   Write distributed boundary-cell-based array.                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type) :: Comm
  integer          :: fh         ! file handle
  real             :: array(:)
  integer          :: disp       ! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: error
!==============================================================================!

  ! Set view for distributed boundary cell data
  ! (this part is the same as in the Read counterpart)
  call Mpi_File_Set_View(fh,                        &
                         disp,                      &
                         MPI_DOUBLE,                &
                         Comm % bnd_cell_map_type,  &
                         'native',                  &
                         MPI_INFO_NULL,             &
                         error)

  ! Write distributed boundary cell data
  call Mpi_File_Write(fh,                 &
                      array,              &
                      Comm % nb_sub,      &
                      MPI_DOUBLE,         &
                      MPI_STATUS_IGNORE,  &
                      error)

  disp = disp + Comm % nb_tot * SIZE_REAL

  end subroutine
