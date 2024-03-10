!==============================================================================!
  subroutine Save_Cell_Real(Backup, Grid, disp, vc, var_name, array)
!------------------------------------------------------------------------------!
!   Writes a vector variable with boundary cells to backup file.               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Backup_Type)      :: Backup
  type(Grid_Type), target :: Grid
  integer(DP)             :: disp
  integer                 :: vc
  character(len=*)        :: var_name
  real                    :: array(-Grid % n_bnd_cells:Grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  type(Comm_Type), pointer :: Comm
  character(SL)            :: vn
  integer                  :: vs  ! variable size
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Backup)
!==============================================================================!

  ! Take alias
  Comm => Grid % Comm

  if(First_Proc()) print *, '# Writing variable: ', trim(var_name)

  ! Increase variable count
  vc = vc + 1

  ! Vector with boundaries
  vn = var_name
  call Comm % Write_Text(fh, vn, disp)

  vs = (Comm % nc_tot + Comm % nb_tot) * RP
  call Comm % Write_Int (fh, vs, disp)

  call Comm % Write_Cell_Real(fh, array(1:Comm % nc_sub), disp)
  call Comm % Write_Bnd_Real (fh, array( -Comm % nb_f:  &
                                         -Comm % nb_l),   disp)

  end subroutine
