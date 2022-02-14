!==============================================================================!
  subroutine Write_Bnd_Real(Comm, fh, array, disp)
!------------------------------------------------------------------------------!
!   Write distributed boundary-cell-based array.                               !
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
                      array,              &
                      Comm % nb_sub,      &
                      comm_type_real,     &
                      MPI_STATUS_IGNORE,  &
                      error)

  disp = disp + Comm % nb_tot * RP

  end subroutine
