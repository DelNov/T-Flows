!==============================================================================!
  subroutine Backup_Mod_Write_Variable(fh, disp, vc, var_name, var)
!------------------------------------------------------------------------------!
!   Writes a whole variable to backup file.                                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer          :: fh, disp, vc
  character(len=*) :: var_name
  type(Var_Type)   :: var
!-----------------------------------[Locals]-----------------------------------!
  type(Comm_Type), pointer :: comm
  character(len=80)        :: vn
  integer                  :: vs  ! variable size
!==============================================================================!

  ! Take alias
  comm => var % pnt_grid % comm

  if(this_proc < 2) print *, '# Writing variable: ', trim(var_name)

  ! Increase variable count
  vc = vc + 1

  ! Vector without boundaries
  vn = var_name
  call Comm_Mod_Write_Text(fh, vn, disp)
  vs = (2*comm % nc_t + 2*comm % nb_t) * SIZE_REAL
  call Comm_Mod_Write_Int (fh, vs, disp)

  call Comm_Mod_Write_Cell_Real(comm, fh, var % n(1:comm % nc_s),   disp)
  call Comm_Mod_Write_Bnd_Real (comm, fh, var % n(-comm % nb_s:-1), disp)
  call Comm_Mod_Write_Bnd_Real (comm, fh, var % q(-comm % nb_s:-1), disp)
  call Comm_Mod_Write_Cell_Real(comm, fh, var % o(1:comm % nc_s),   disp)

  end subroutine
