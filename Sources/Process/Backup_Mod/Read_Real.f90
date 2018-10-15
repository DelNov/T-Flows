!==============================================================================!
  subroutine Backup_Mod_Read_Real(fh, disp, vc, var_name, var_value)
!------------------------------------------------------------------------------!
!   Reads a single named real variable from backup file.                      !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Comm_Mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer          :: fh, disp, vc
  character(len=*) :: var_name
  real             :: var_value
!-----------------------------------[Locals]-----------------------------------!
  character(len=80) :: vn
  integer           :: vo, disp_loop, cnt_loop
!==============================================================================!

  cnt_loop  = 0
  disp_loop = 0

  !--------------------------------------------------------!
  !   Browse the entire file until you find the variable   !
  !--------------------------------------------------------!
  do

    call Comm_Mod_Read_Text(fh, vn, disp_loop)  ! variable name
    call Comm_Mod_Read_Int (fh, vo, disp_loop)  ! variable offset

    ! If variable is found, read it and retrun
    if(vn .eq. var_name) then
      if(this_proc < 2) print *, '# Reading variable: ', trim(vn)
      call Comm_Mod_Read_Real(fh, var_value, disp_loop)
      disp = disp_loop
      return

    ! If variable not found, advance the offset only
    else
      disp_loop = disp_loop + vo
    end if

    ! Check if variable is in the file
    cnt_loop = cnt_loop + 1
    if(cnt_loop > vc) goto 1

  end do

1 if(this_proc < 2) print *, '# Variable: ', trim(vn), ' not found!'

  end subroutine
