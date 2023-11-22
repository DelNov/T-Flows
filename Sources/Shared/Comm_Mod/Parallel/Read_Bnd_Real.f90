!==============================================================================!
  subroutine Read_Bnd_Real(Comm, fh, array, disp)
!------------------------------------------------------------------------------!
!   Read distributed boundary-cell-based array.                                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type), intent(in)    :: Comm
  type(Mpi_File),   intent(in)    :: fh         ! file handle
  real,             intent(out)   :: array(:)
  integer(DP),      intent(inout) :: disp       ! displacement in bytes
!-----------------------------------[Locals]-----------------------------------!
  integer :: error = 0
!==============================================================================!

  ! Set view for distributed boundary cell data
  ! (this part is the same as in Write counterpart)
  call Mpi_File_Set_View(fh,                        &
                         disp,                      &
                         comm_type_real,            &
                         Comm % bnd_cell_map_type,  &
                         'native',                  &
                         MPI_INFO_NULL,             &
                         error)

  ! Read distributed boundary cell data 
  call Mpi_File_Read(fh,                 &
                     array,              &
                     Comm % nb_sub,      &
                     comm_type_real,     &
                     MPI_STATUS_IGNORE,  &
                     error)

  disp = disp + Comm % nb_tot * RP

  end subroutine
