!==============================================================================!
  subroutine Save_Log(Backup, Comm, disp, vc, var_name, var_value)
!------------------------------------------------------------------------------!
!   Writes a named logical variable to backup file.                            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Backup_Type) :: Backup
  type(Comm_Type)    :: Comm
  integer(DP)        :: disp
  integer            :: vc
  character(len=*)   :: var_name
  logical            :: var_value
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: vn
  integer       :: vs  ! variable size
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Backup)
!==============================================================================!

  if(First_Proc()) print *, '# Writing variable: ', trim(var_name)

  ! Increase variable count
  vc = vc + 1

  ! Just store one named logical
  vn = var_name;  call Comm % Write_Text(fh, vn, disp)
  vs = LP;  call Comm % Write_Int (fh, vs, disp)

  call Comm % Write_Log(fh, var_value, disp)

  end subroutine
