!==============================================================================!
  subroutine Write_Face_Real(Comm, fh, array, disp)
!------------------------------------------------------------------------------!
!   Write distributed face-based array.                                        !
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

  ! Set view for distributed face data
  ! (this part is the same as in Read counterpart)
  call Mpi_File_Set_View(fh,                    &
                         disp,                  &
                         MPI_DOUBLE,            &
                         Comm % face_map_type,  &
                         'native',              &
                         MPI_INFO_NULL,         &
                         error)

  ! Write distributed face data
  call Mpi_File_Write(fh,                 &
                      array,              &
                      Comm % nf_sub,      &
                      MPI_DOUBLE,         &
                      MPI_STATUS_IGNORE,  &
                      error)

  disp = disp + Comm % nf_tot * SIZE_REAL

  end subroutine
