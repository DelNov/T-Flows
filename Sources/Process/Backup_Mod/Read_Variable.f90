!==============================================================================!
  subroutine Backup_Mod_Read_Variable(fh, disp, var_name, var)
!------------------------------------------------------------------------------!
!   Reads a whole variable from backup file.                                   !
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
  type(Var_Type)   :: var
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
    if(vn .eq. var_name) then
      if(this_proc < 2) print *, '# Reading variable: ', trim(vn)
      call Comm_Mod_Read_Cell_Real(fh, var % n(1:nc_s),   disp_loop)
      call Comm_Mod_Read_Bnd_Real (fh, var % n(-nb_s:-1), disp_loop)
      call Comm_Mod_Read_Bnd_Real (fh, var % q(-nb_s:-1), disp_loop)
      call Comm_Mod_Read_Cell_Real(fh, var % o  (1:nc_s), disp_loop)
      call Comm_Mod_Read_Cell_Real(fh, var % a  (1:nc_s), disp_loop)
      call Comm_Mod_Read_Cell_Real(fh, var % a_o(1:nc_s), disp_loop)
      call Comm_Mod_Read_Cell_Real(fh, var % c  (1:nc_s), disp_loop)
      call Comm_Mod_Read_Cell_Real(fh, var % c_o(1:nc_s), disp_loop)
      call Comm_Mod_Read_Cell_Real(fh, var % d_o(1:nc_s), disp_loop)
      disp = disp_loop
      return

    ! If variable not found, advance the offset only
    else
      disp_loop = disp_loop + vs
    end if

  end do

  if(this_proc < 2) print *, '# Variable: ', trim(vn), ' not found!'

  end subroutine
