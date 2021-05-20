!==============================================================================!
  subroutine Backup_Mod_Write_Cell_Real(Grid, fh, disp, vc, var_name, array)
!------------------------------------------------------------------------------!
!   Writes a vector variable with boundary cells to backup file.               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type), target :: Grid
  integer                 :: fh, disp, vc
  character(len=*)        :: var_name
  real                    :: array(-Grid % n_bnd_cells:Grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  type(Comm_Type), pointer :: Comm
  character(SL)            :: vn
  integer                  :: vs  ! variable size
!==============================================================================!

  ! Take alias
  Comm => Grid % Comm

  if(this_proc < 2) print *, '# Writing variable: ', trim(var_name)

  ! Increase variable count
  vc = vc + 1

  ! Vector with boundaries
  vn = var_name
  call Comm % Write_Text(fh, vn, disp)

  vs = (Comm % nc_tot + Comm % nb_tot) * SIZE_REAL
  call Comm % Write_Int (fh, vs, disp)

  call Comm % Write_Cell_Real(fh, array(1:Comm % nc_sub), disp)
  call Comm % Write_Bnd_Real (fh, array( -Comm % nb_f:  &
                                         -Comm % nb_l),   disp)

  end subroutine
