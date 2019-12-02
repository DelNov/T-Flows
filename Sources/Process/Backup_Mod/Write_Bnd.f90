!==============================================================================!
  subroutine Backup_Mod_Write_Bnd(comm, fh, disp, vc, var_name, array)
!------------------------------------------------------------------------------!
!   Writes a vector variable with boundary cells to backup file.               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Comm_Type)  :: comm
  integer          :: fh, disp, vc
  character(len=*) :: var_name
  real             :: array(-comm % nb_s:-1)
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
  vs = comm % nb_t * SIZE_REAL
  call Comm_Mod_Write_Int (fh, vs, disp)

  call Comm_Mod_Write_Bnd_Real(comm, fh, array(-comm % nb_s:-1), disp)

  end subroutine
