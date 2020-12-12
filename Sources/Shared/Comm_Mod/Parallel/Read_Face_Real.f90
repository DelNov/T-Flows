!==============================================================================!
  subroutine Comm_Mod_Read_Face_Real(comm, fh, array, disp)
!------------------------------------------------------------------------------!
!   Read distributed face-based array.                                         !
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

  ! Set view for distributed face data
  ! (this part is the same as in Write counterpart)
  call Mpi_File_Set_View(fh,             &
                         disp,           &
                         MPI_DOUBLE,     &
                         face_map_type,  &
                         'native',       &
                         MPI_INFO_NULL,  &
                         error)

  ! Read distributed face data
  call Mpi_File_Read(fh,                 &
                     array,              &
                     comm % nf_sub,      &
                     MPI_DOUBLE,         &
                     MPI_STATUS_IGNORE,  &
                     error)

  disp = disp + comm % nf_tot * SIZE_REAL

  end subroutine
