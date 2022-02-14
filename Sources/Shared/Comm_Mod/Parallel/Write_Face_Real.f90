!==============================================================================!
  subroutine Write_Face_Real(Comm, fh, array, disp)
!------------------------------------------------------------------------------!
!   Write distributed face-based array.  (Not used any more.)                  !
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

  ! Set view for distributed face data
  ! (this part is the same as in Read counterpart)
  call Mpi_File_Set_View(fh,                    &
                         disp,                  &
                         comm_type_real,        &
                         Comm % face_map_type,  &
                         'native',              &
                         MPI_INFO_NULL,         &
                         error)

  ! Write distributed face data
  call Mpi_File_Write(fh,                 &
                      array,              &
                      Comm % nf_sub,      &
                      comm_type_real,     &
                      MPI_STATUS_IGNORE,  &
                      error)

  disp = disp + Comm % nf_tot * RP

  end subroutine
