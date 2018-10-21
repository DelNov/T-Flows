!==============================================================================!
  subroutine Backup_Mod_Write_Variable(fh, disp, vc, var_name, var)
!------------------------------------------------------------------------------!
!   Writes a whole variable to backup file.                                    !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Comm_Mod
  use Grid_Mod
  use Var_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer          :: fh, disp, vc
  character(len=*) :: var_name
  type(Var_Type)   :: var
!-----------------------------------[Locals]-----------------------------------!
  character(len=80) :: vn
  integer           :: vs  ! variable size
!==============================================================================!

  if(this_proc < 2) print *, '# Writing variable: ', trim(var_name)

  ! Increase variable count
  vc = vc + 1

  ! Vector without boundaries
  vn = var_name;                      call Comm_Mod_Write_Text(fh, vn, disp)
  vs = (2*nc_t + 2*nb_t) * SIZE_REAL; call Comm_Mod_Write_Int (fh, vs, disp)

  call Comm_Mod_Write_Cell_Real(fh, var % n(1:nc_s),   disp)
  call Comm_Mod_Write_Bnd_Real (fh, var % n(-nb_s:-1), disp)
  call Comm_Mod_Write_Bnd_Real (fh, var % q(-nb_s:-1), disp)
  call Comm_Mod_Write_Cell_Real(fh, var % o(1:nc_s),   disp)

  end subroutine
