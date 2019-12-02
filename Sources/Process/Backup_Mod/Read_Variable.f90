!==============================================================================!
  subroutine Backup_Mod_Read_Variable(fh, disp, vc, var_name, var)
!------------------------------------------------------------------------------!
!   Reads a whole variable from backup file.                                   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer          :: fh, disp, vc
  character(len=*) :: var_name
  type(Var_Type)   :: var
!-----------------------------------[Locals]-----------------------------------!
  type(Comm_Type), pointer :: comm
  character(len=80)        :: vn
  integer                  :: vs, disp_loop, cnt_loop
!==============================================================================!

  ! Take alias
  comm => var % pnt_grid % comm

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
      call Comm_Mod_Read_Cell_Real(comm,fh,var % n(1:comm % nc_s),  disp_loop)
      call Comm_Mod_Read_Bnd_Real (comm,fh,var % n(-comm % nb_s:-1),disp_loop)
      call Comm_Mod_Read_Bnd_Real (comm,fh,var % q(-comm % nb_s:-1),disp_loop)
      call Comm_Mod_Read_Cell_Real(comm,fh,var % o(1:comm % nc_s),  disp_loop)
      disp = disp_loop
      return

    ! If variable not found, advance the offset only
    else
      disp_loop = disp_loop + vs
    end if

    ! Check if variable is in the file
    if(cnt_loop >= vc) goto 1

  end do

1 if(this_proc < 2) print *, '# Variable: ', trim(var_name), ' not found!',  &
                             'Setting the values to 0.0!'
  var % n(:) = 0.0
  var % q(:) = 0.0
  var % o(:) = 0.0

  end subroutine
