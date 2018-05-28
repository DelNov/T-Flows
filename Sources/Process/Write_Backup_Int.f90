!==============================================================================!
  subroutine Write_Backup_Int(fh, disp, var_name, var_value)
!------------------------------------------------------------------------------!
!   Writes a single named integer variable to backup file.                     !
!----------------------------------[Modules]-----------------------------------!
  use Comm_Mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer          :: fh, disp
  character(len=*) :: var_name
  integer          :: var_value
!-----------------------------------[Locals]-----------------------------------!
  character(len=80) :: vn
  integer           :: vs  ! variable size
!==============================================================================!

  if(this_proc < 2) print *, '# Writing variable: ', trim(var_name)

  ! Just store one named integer
  vn = var_name;  call Comm_Mod_Write_Text(fh, vn, disp)
  vs = SIZE_INT;  call Comm_Mod_Write_Int (fh, vs, disp)

  call Comm_Mod_Write_Int (fh, var_value, disp)

  end subroutine
