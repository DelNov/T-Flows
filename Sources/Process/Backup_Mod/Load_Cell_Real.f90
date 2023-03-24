!==============================================================================!
  subroutine Load_Cell_Real(Backup, Grid, disp, vc, var_name, array)
!------------------------------------------------------------------------------!
!   Reads a vector variable with boundary cells from a backup file.            !
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
  integer                  :: vs, cnt_loop, nb, nc
  integer(DP)              :: disp_loop
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Backup)
!==============================================================================!

  ! Take alias
  Comm => Grid % Comm
  nb = Grid % n_bnd_cells
  nc = Grid % n_cells

  cnt_loop  = 0
  disp_loop = 0

  !--------------------------------------------------------!
  !   Browse the entire file until you find the variable   !
  !--------------------------------------------------------!
  do

    ! Increase counter
    cnt_loop = cnt_loop + 1

    call Comm % Read_Text(fh, vn, disp_loop)  ! variable name
    call Comm % Read_Int (fh, vs, disp_loop)  ! variable size

    ! If variable is found, read it and retrun
    if(vn .eq. var_name) then
      if(First_Proc()) print *, '# Reading variable: ', trim(vn)
      call Comm % Read_Cell_Real(fh, array(1:Comm % nc_sub), disp_loop)
      call Comm % Read_Bnd_Real (fh, array( -Comm % nb_f:  &
                                            -Comm % nb_l),   disp_loop)
      call Grid % Exchange_Cells_Real(array(-nb:nc))
      disp = disp_loop
      return

    ! If variable not found, advance the offset only
    else
      disp_loop = disp_loop + vs
    end if

    ! Check if variable is in the file
    if(cnt_loop >= vc) goto 1

  end do

1 if(First_Proc()) print *, '# Variable: ', trim(var_name), ' not found! ',  &
                             'Setting the values to 0.0!'
  array(:) = 0.0

  end subroutine
