!==============================================================================!
  subroutine Load_Variable(Backup, disp, vc, var_name, Flow, var)
!------------------------------------------------------------------------------!
!   Reads a whole variable from backup file.                                   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Backup_Type) :: Backup
  type(Field_Type)   :: Flow
  integer(DP)        :: disp
  integer            :: vc
  character(len=*)   :: var_name
  type(Var_Type)     :: var
!-----------------------------------[Locals]-----------------------------------!
  type(Comm_Type), pointer :: Comm
  type(Grid_Type), pointer :: Grid
  character(SL)            :: vn
  integer                  :: vs, cnt_loop, nb, nc
  integer(DP)              :: disp_loop
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Backup)
!==============================================================================!

  ! Take aliases
  Grid => var % pnt_grid
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
      call Comm % Read_Cell_Real(fh,var % n(1:Comm % nc_sub), disp_loop)
      call Comm % Read_Bnd_Real (fh,var % n( -Comm % nb_f:  &
                                             -Comm % nb_l),   disp_loop)
      call Comm % Read_Cell_Real(fh,var % q(1:Comm % nc_sub), disp_loop)
      call Comm % Read_Bnd_Real (fh,var % q( -Comm % nb_f:  &
                                             -Comm % nb_l),   disp_loop)
      call Comm % Read_Cell_Real(fh,var % o(1:Comm % nc_sub), disp_loop)
      call Comm % Read_Bnd_Real (fh,var % o( -Comm % nb_f:  &
                                             -Comm % nb_l),   disp_loop)

      ! Refresh buffers for "n", "q" and "o" (not sure if needed)
      call Grid % Exchange_Cells_Real(var % n(-nb:nc))
      call Grid % Exchange_Cells_Real(var % q(-nb:nc))
      call Grid % Exchange_Cells_Real(var % o(-nb:nc))

      ! Compute fresh gradients
      call Flow % Grad_Variable(var)

      disp = disp_loop
      return

    ! If variable not found, advance the offset only
    else
      disp_loop = disp_loop + vs
    end if

    ! Check if variable is in the file
    if(cnt_loop >= vc) goto 1

  end do

1 if(First_Proc()) print *, '# Variable: ', trim(var_name), ' not found!',  &
                             'Setting the values to 0.0!'
  var % n(:) = 0.0
  var % q(:) = 0.0
  var % o(:) = 0.0

  end subroutine
