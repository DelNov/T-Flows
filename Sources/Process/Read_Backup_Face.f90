!==============================================================================!
  subroutine Read_Backup_Face(fh, disp, var_name, fcom1)
!------------------------------------------------------------------------------!
!   Reads a vector variable with boundary cells from a backup file.            !
!------------------------------------------------------------------------------!
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
  integer           :: vs, disp_loop
!==============================================================================!

  disp_loop = 0

  !--------------------------------------------------------!
  !   Browse the entire file until you find the variable   !
  !--------------------------------------------------------!
  do

    call Comm_Mod_Read_Text(fh, vn, disp_loop)  ! variable name
    call Comm_Mod_Read_Int (fh, vs, disp_loop)  ! variable size  

    ! If variable is found, read it and retrun
    if(vn .eq. var_name) then
      if(this_proc < 2) print *, '# Reading variable: ', trim(vn)
      call Comm_Mod_Read_Face_Real(fh, fcom1(1:nf_s+nbf_s), disp_loop)
      disp = disp_loop

      ! Fix signs for the face fcom1es back
      do s = nf_s + 1, nf_s + nbf_s
        fcom1(s) = fcom1(s) * buf_face_sgn(s-nf_s)
      end do

      return

    ! If variable not found, advance the offset only
    else
      disp_loop = disp_loop + vs
    end if

  end do

  if(this_proc < 2) print *, '# Variable: ', trim(vn), ' not found!'

  end subroutine
