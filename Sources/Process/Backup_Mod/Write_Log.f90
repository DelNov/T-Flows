!==============================================================================!
  subroutine Backup_Mod_Write_Log(Comm, fh, disp, vc, var_name, var_value)
!------------------------------------------------------------------------------!
!   Writes a named logical variable to backup file.                            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Comm_Type)  :: Comm
  integer          :: fh, disp, vc
  character(len=*) :: var_name
  logical          :: var_value
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: vn
  integer       :: vs  ! variable size
!==============================================================================!

  if(this_proc < 2) print *, '# Writing variable: ', trim(var_name)

  ! Increase variable count
  vc = vc + 1

  ! Just store one named logical
  vn = var_name;  call Comm % Write_Text(fh, vn, disp)
  vs = SIZE_LOG;  call Comm % Write_Int (fh, vs, disp)

  call Comm % Write_Log(fh, var_value, disp)

  end subroutine
