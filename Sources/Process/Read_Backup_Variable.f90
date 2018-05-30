!==============================================================================!
  subroutine Read_Backup_Variable(fh, disp, var_name, var1)
!------------------------------------------------------------------------------!
!   Reads a vector variable with boundary cells from a backup file.            !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Comm_Mod
  use Grid_Mod
  use Var_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer          :: fh, disp
  character(len=*) :: var_name
  type(Var_Type)   :: var1
!-----------------------------------[Locals]-----------------------------------!
  character(len=80) :: vn
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
    if(vn == var_name) then
      if(this_proc < 2) print *, '# Reading variable: ', trim(vn)
      call Comm_Mod_Read_Cell_Real(fh, var1 % n(1:nc_s),   disp_loop)
      call Comm_Mod_Read_Bnd_Real (fh, var1 % n(-nb_s:-1), disp_loop)
      call Comm_Mod_Read_Bnd_Real (fh, var1 % q(-nb_s:-1), disp_loop)
      call Comm_Mod_Read_Cell_Real(fh, var1 % o  (1:nc_s), disp_loop)
      call Comm_Mod_Read_Cell_Real(fh, var1 % a  (1:nc_s), disp_loop)
      call Comm_Mod_Read_Cell_Real(fh, var1 % a_o(1:nc_s), disp_loop)
      call Comm_Mod_Read_Cell_Real(fh, var1 % c  (1:nc_s), disp_loop)
      call Comm_Mod_Read_Cell_Real(fh, var1 % c_o(1:nc_s), disp_loop)
      call Comm_Mod_Read_Cell_Real(fh, var1 % d_o(1:nc_s), disp_loop)
      disp = disp_loop
      return

    ! If variable not found, advance the offset only
    else
      disp_loop = disp_loop + vs
    end if

  end do

  if(this_proc < 2) print *, '# Variable: ', trim(vn), ' not found!'

  end subroutine
