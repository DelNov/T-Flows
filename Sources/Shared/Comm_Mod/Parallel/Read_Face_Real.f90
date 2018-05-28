!==============================================================================!
  subroutine Comm_Mod_Read_Face_Real(fh, array, disp)
!------------------------------------------------------------------------------!
!   Read distributed face-based array.                                         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: fh         ! file handle
  real    :: array(:)
  integer :: disp       ! displacement
!-----------------------------------[Locals]-----------------------------------!
  integer :: s, error
!==============================================================================!

  ! Read arrays only for parallel runs
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

    ! Read distributed face data 
    call Mpi_File_Read(fh,                 &
                       buf_face_val,       &
                       nbf_s,              &   
                       MPI_DOUBLE,         &   
                       MPI_STATUS_IGNORE,  &
                       error)
    do s = 1, nbf_s
      array(nf_s + s) = buf_face_val(buf_face_ord(s))
    end do

  end if

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
                     face_val,           &
                     nf_s,               &   
                     MPI_DOUBLE,         &   
                     MPI_STATUS_IGNORE,  &
                     error)
  do s = 1, nf_s
    array(s) = face_val(face_ord(s))
  end do

  disp = disp + nf_t * SIZE_REAL

  end subroutine
