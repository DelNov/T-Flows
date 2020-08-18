!==============================================================================!
  subroutine Backup_Mod_Write_Real(fh, disp, vc, var_name, var_value)
!------------------------------------------------------------------------------!
!   Writes a single named real variable to backup file.                        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer          :: fh, disp, vc
  character(len=*) :: var_name
  real             :: var_value
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: vn
  integer       :: vs  ! variable size
!==============================================================================!

  if(this_proc < 2) print *, '# Writing variable: ', trim(var_name)

  ! Increase variable count
  vc = vc + 1

  ! Just store one named real number
  vn = var_name;  call Comm_Mod_Write_Text(fh, vn, disp)
  vs = SIZE_REAL; call Comm_Mod_Write_Int (fh, vs, disp)

  call Comm_Mod_Write_Real(fh, var_value, disp)

  end subroutine
