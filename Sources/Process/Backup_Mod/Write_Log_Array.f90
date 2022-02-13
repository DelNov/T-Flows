!==============================================================================!
  subroutine Backup_Mod_Write_Log_Array(Comm, disp, vc, arr_name, arr_value)
!------------------------------------------------------------------------------!
!   Writes a named logical array to backup file.                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Comm_Type)       :: Comm
  integer(DP)           :: disp
  integer               :: vc
  character(len=*)      :: arr_name
  logical, dimension(:) :: arr_value
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: vn
  integer       :: vs, length  ! variable size
!==============================================================================!

  length = size(arr_value)

  if(this_proc < 2) print *, '# Writing array: ', trim(arr_name)

  ! Increase variable count
  vc = vc + 1

  ! Just store one named integer array
  vn = arr_name;           call Comm % Write_Text(fh, vn, disp)
  vs = length * SIZE_LOG;  call Comm % Write_Int (fh, vs, disp)

  call Comm % Write_Log_Array(fh, arr_value(1:length), disp)

  end subroutine
