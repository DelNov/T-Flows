!==============================================================================!
  subroutine Backup_Mod_Write_Int(fh, disp, vc, var_name, var_value)
!------------------------------------------------------------------------------!
!   Writes a single named integer variable to backup file.                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer          :: fh, disp, vc
  character(len=*) :: var_name
  integer          :: var_value
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: vn
  integer       :: vs  ! variable size
!==============================================================================!

  if(this_proc < 2) print *, '# Writing variable: ', trim(var_name)

  ! Increase variable count
  vc = vc + 1

  ! Just store one named integer
  vn = var_name;  call Comm_Mod_Write_Text(fh, vn, disp)
  vs = SIZE_INT;  call Comm_Mod_Write_Int (fh, vs, disp)

  call Comm_Mod_Write_Int(fh, var_value, disp)

  end subroutine
