!==============================================================================!
  subroutine Backup_Mod_Write_Cell_Real(grid, fh, disp, vc, var_name, array)
!------------------------------------------------------------------------------!
!   Writes a vector variable with boundary cells to backup file.               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type), target :: grid
  integer                 :: fh, disp, vc
  character(len=*)        :: var_name
  real                    :: array(-grid % n_bnd_cells:grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  type(Comm_Type), pointer :: comm
  character(SL)            :: vn
  integer                  :: vs  ! variable size
!==============================================================================!

  ! Take alias
  comm => grid % comm

  if(this_proc < 2) print *, '# Writing variable: ', trim(var_name)

  ! Increase variable count
  vc = vc + 1

  ! Vector with boundaries
  vn = var_name
  call Comm_Mod_Write_Text(fh, vn, disp)

  vs = (comm % nc_t + comm % nb_t) * SIZE_REAL
  call Comm_Mod_Write_Int (fh, vs, disp)

  call Comm_Mod_Write_Cell_Real(comm, fh, array(1:comm % nc_s),   disp)
  call Comm_Mod_Write_Bnd_Real (comm, fh, array(-comm % nb_f:  &
                                                -comm % nb_l), disp)

  end subroutine
