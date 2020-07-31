!==============================================================================!
  subroutine Backup_Mod_Read_Cell_Real(grid, fh, disp, vc, var_name, array)
!------------------------------------------------------------------------------!
!   Reads a vector variable with boundary cells from a backup file.            !
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
  integer                  :: vs, disp_loop, cnt_loop, nb, nc
!==============================================================================!

  ! Take alias
  comm => grid % comm
  nb = grid % n_bnd_cells
  nc = grid % n_cells

  cnt_loop  = 0
  disp_loop = 0

  !--------------------------------------------------------!
  !   Browse the entire file until you find the variable   !
  !--------------------------------------------------------!
  do

    ! Increase counter
    cnt_loop = cnt_loop + 1

    call Comm_Mod_Read_Text(fh, vn, disp_loop)  ! variable name
    call Comm_Mod_Read_Int (fh, vs, disp_loop)  ! variable size  

    ! If variable is found, read it and retrun
    if(vn .eq. var_name) then
      if(this_proc < 2) print *, '# Reading variable: ', trim(vn)
      call Comm_Mod_Read_Cell_Real(comm, fh, array(1:comm % nc_s),   disp_loop)
      call Comm_Mod_Read_Bnd_Real (comm, fh, array(-comm % nb_f:  &
                                                   -comm % nb_l), disp_loop)
      call Grid_Mod_Exchange_Cells_Real(grid, array(-nb:nc))
      disp = disp_loop
      return

    ! If variable not found, advance the offset only
    else
      disp_loop = disp_loop + vs
    end if

    ! Check if variable is in the file
    if(cnt_loop >= vc) goto 1

  end do

1 if(this_proc < 2) print *, '# Variable: ', trim(var_name), ' not found! ',  &
                             'Setting the values to 0.0!'
  array(:) = 0.0

  end subroutine
