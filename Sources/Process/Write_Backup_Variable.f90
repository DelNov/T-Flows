!==============================================================================!
  subroutine Write_Backup_Variable(fh, disp, var_name, var1)
!------------------------------------------------------------------------------!
!   Writes a whole variable cells to backup file.                              !
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
  integer           :: vs  ! variable size
!==============================================================================!

  if(this_proc < 2) print *, '# Writing variable: ', trim(var_name)

  ! Vector without boundaries
  vn = var_name;                      call Comm_Mod_Write_Text(fh, vn, disp)
  vs = (7*nc_t + 2*nb_t) * SIZE_REAL; call Comm_Mod_Write_Int (fh, vs, disp)

  call Comm_Mod_Write_Cell_Real(fh, var1 % n(1:nc_s),   disp)
  call Comm_Mod_Write_Bnd_Real (fh, var1 % n(-nb_s:-1), disp)
  call Comm_Mod_Write_Bnd_Real (fh, var1 % q(-nb_s:-1), disp)
  call Comm_Mod_Write_Cell_Real(fh, var1 % o  (1:nc_s), disp)
  call Comm_Mod_Write_Cell_Real(fh, var1 % a  (1:nc_s), disp)
  call Comm_Mod_Write_Cell_Real(fh, var1 % a_o(1:nc_s), disp)
  call Comm_Mod_Write_Cell_Real(fh, var1 % c  (1:nc_s), disp)
  call Comm_Mod_Write_Cell_Real(fh, var1 % c_o(1:nc_s), disp)
  call Comm_Mod_Write_Cell_Real(fh, var1 % d_o(1:nc_s), disp)

  end subroutine
