!==============================================================================!
  subroutine Backup_Mod_Read_Face_Real(Grid, fh, disp, vc, var_name, array,  &
                                       correct_sign)
!------------------------------------------------------------------------------!
!   Reads a vector variable with face values from a backup file.               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type), target :: Grid
  integer                 :: fh, disp, vc
  character(len=*)        :: var_name
  real                    :: array(Grid % n_faces)
  logical, optional       :: correct_sign  ! for face fluxes, signs might have
                                           ! to be changed (check it one day)
!-----------------------------------[Locals]-----------------------------------!
  type(Comm_Type), pointer :: Comm
  character(len=80)        :: vn
  integer                  :: vs, disp_loop, cnt_loop, nf
  integer                  :: s, c1, c2, cg1, cg2
!==============================================================================!

  ! Take alias
  Comm => Grid % Comm
  nf = Grid % n_faces

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
      if(this_proc < 2) print *, '# Reading variable: ', trim(vn)
      call Comm % Read_Face_Real(fh, array(1:Comm % nf_sub), disp_loop)
      disp = disp_loop

      ! Correct the signs of fluxes  (Remember, they are defined
      ! defined to be positive from cg1 to cg2; and cg2 > cg1)
      if(present(correct_sign)) then
        if(correct_sign) then
          do s = 1, Grid % n_faces
            c1  = Grid % faces_c(1,s)
            c2  = Grid % faces_c(2,s)
            cg1 = Grid % Comm % cell_glo(c1)
            cg2 = Grid % Comm % cell_glo(c2)
            if(cg2 > 0 .and. cg2 < cg1) then
              array(s) = -array(s)
            end if
          end do
        end if
      end if

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
