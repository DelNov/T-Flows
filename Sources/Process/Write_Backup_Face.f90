!==============================================================================!
  subroutine Write_Backup_Face(fh, disp, var_name, fcom1)
!------------------------------------------------------------------------------!
!   Writes a vector variable with boundary cells to backup file.               !
!----------------------------------[Modules]-----------------------------------!
  use Comm_Mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer          :: fh, disp
  character(len=*) :: var_name
  real             :: fcom1(1:nf_s+nbf_s)
!-----------------------------------[Locals]-----------------------------------!
  character(len=80) :: vn
  integer           :: s
  integer           :: vs  ! variable size
!==============================================================================!

  if(this_proc < 2) print *, '# Writing variable: ', trim(var_name)

  ! Change the sign of fcom1 where necessary
  do s = nf_s + 1, nf_s + nbf_s
    fcom1(s) = fcom1(s) * buf_face_sgn(s-nf_s)
  end do

  ! Perform the actual saving
  vn = var_name;         call Comm_Mod_Write_Text(fh, vn, disp)
  vs = nf_t * SIZE_REAL; call Comm_Mod_Write_Int (fh, vs, disp)

  call Comm_Mod_Write_Face_Real(fh, fcom1(1:nf_s+nbf_s), disp)

  ! Fix signs for the face fcom1 back
  do s = nf_s + 1, nf_s + nbf_s
    fcom1(s) = fcom1(s) * buf_face_sgn(s-nf_s)
  end do


  end subroutine
