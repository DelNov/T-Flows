!==============================================================================!
  subroutine Backup_Mod_Write_Variable_Mean(fh, disp, var_name, var)
!------------------------------------------------------------------------------!
!   Writes variable's mean with boundary cells from a backup file.             !
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
  integer           :: vs  ! variable size
!==============================================================================!

  if(this_proc < 2) print *, '# Writing variable: ', trim(var_name)

  ! Vector without boundaries
  vn = var_name;                  call Comm_Mod_Write_Text(fh, vn, disp)
  vs = (nc_t + nb_t) * SIZE_REAL; call Comm_Mod_Write_Int (fh, vs, disp)

  call Comm_Mod_Write_Cell_Real(fh, var % mean(1:nc_s),   disp)
  call Comm_Mod_Write_Bnd_Real (fh, var % mean(-nb_s:-1), disp)

  end subroutine
