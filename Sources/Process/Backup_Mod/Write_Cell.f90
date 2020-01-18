!==============================================================================!
  subroutine Backup_Mod_Write_Cell(comm, fh, disp, vc, var_name, array)
!------------------------------------------------------------------------------!
!   Writes a vector variable with boundary cells to backup file.               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Comm_Type)  :: comm
  integer          :: fh, disp, vc
  character(len=*) :: var_name
  real             :: array(1:comm % nc_s)
!-----------------------------------[Locals]-----------------------------------!
  character(len=80) :: vn
  integer           :: vs  ! variable size
!==============================================================================!

  if(this_proc < 2) print *, '# Writing variable: ', trim(var_name)

  ! Increase variable count
  vc = vc + 1

  ! Vector with boundaries
  vn = var_name
  call Comm_Mod_Write_Text(fh, vn, disp)
  vs = comm % nc_t * SIZE_REAL
  call Comm_Mod_Write_Int (fh, vs, disp)

  call Comm_Mod_Write_Cell_Real(comm, fh, array(1:comm % nc_s), disp)

  end subroutine
