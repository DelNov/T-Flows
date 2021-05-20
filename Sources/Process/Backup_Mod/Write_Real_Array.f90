!==============================================================================!
  subroutine Backup_Mod_Write_Real_Array(Comm, fh, disp,  &
                                         vc, arr_name, arr_value)
!------------------------------------------------------------------------------!
!   Writes a name real array to backup file.                                   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Comm_Type)    :: Comm
  integer            :: fh, disp, vc
  character(len=*)   :: arr_name
  real, dimension(:) :: arr_value
!-----------------------------------[Locals]-----------------------------------!
  character(SL)   :: vn
  integer         :: vs, length  ! variable size
!==============================================================================!

  length = size(arr_value)

  if(this_proc < 2) print *, '# Writing array: ', trim(arr_name)

  ! Increase variable count
  vc = vc + 1

  ! Just store one named real number
  vn = arr_name;           call Comm % Write_Text(fh, vn, disp)
  vs = length * SIZE_REAL; call Comm % Write_Int (fh, vs, disp)

  call Comm % Write_Real_Array(fh, arr_value(1:length), disp)

  end subroutine
