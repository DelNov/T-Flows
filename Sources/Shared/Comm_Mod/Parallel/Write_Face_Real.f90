!==============================================================================!
  subroutine Comm_Mod_Write_Face_Real(comm, fh, array, disp)
!------------------------------------------------------------------------------!
!   Write distributed face-based array.                                        !
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
  ! (this part is the same as in Read counterpart)
  call Mpi_File_Set_View(fh,             &
                         disp,           &
                         MPI_DOUBLE,     &
                         face_map_type,  &
                         'native',       &
                         MPI_INFO_NULL,  &
                         error)

  ! Write distributed face data
  call Mpi_File_Write(fh,                 &
                      array,              &
                      comm % nf_s,        &
                      MPI_DOUBLE,         &
                      MPI_STATUS_IGNORE,  &
                      error)

  disp = disp + comm % nf_t * SIZE_REAL

  end subroutine
