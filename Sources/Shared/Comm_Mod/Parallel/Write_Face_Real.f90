!==============================================================================!
  subroutine Comm_Mod_Write_Face_Real(fh, array, disp)
!------------------------------------------------------------------------------!
!   Write distributed face-based array.                                        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: fh         ! file handle
  real    :: array(:)
  integer :: disp       ! displacement
!-----------------------------------[Locals]-----------------------------------!
  integer :: s, error
!==============================================================================!

  ! Save arrays only for parallel runs
  if(n_proc > 1) then

    ! Set view for distributed array face data 
    ! (this part is the same as in Read counterpart)
    call Mpi_File_Set_View(fh,                 &   
                           disp,               &   
                           MPI_DOUBLE,         &
                           buf_face_map_type,  &   
                           'native',           &   
                           MPI_INFO_NULL,      &
                           error)

    ! Write distributed face data 
    do s = 1, nbf_s
      buf_face_val(s) = array( nf_s + buf_face_ord(s) )
    end do
    call Mpi_File_Write(fh,                 &
                        buf_face_val,       &
                        nbf_s,              &   
                        MPI_DOUBLE,         &   
                        MPI_STATUS_IGNORE,  &
                        error)

  end if

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
  do s = 1, nf_s
    face_val(s) = array( face_ord(s) )
  end do
  call Mpi_File_Write(fh,                 &
                      face_val,           &
                      nf_s,               &   
                      MPI_DOUBLE,         &   
                      MPI_STATUS_IGNORE,  &
                      error)

  disp = disp + nf_t * SIZE_REAL

  end subroutine
