!==============================================================================!
  subroutine Load_Face_Real(Backup, Grid, disp, vc, var_name, array, corr_sign)
!------------------------------------------------------------------------------!
!   Reads a vector variable with face values from a backup file.               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Backup_Type)      :: Backup
  type(Grid_Type), target :: Grid
  integer                 :: disp, vc
  character(len=*)        :: var_name
  real                    :: array(Grid % n_faces)
  logical, optional       :: corr_sign  ! for face fluxes, signs might have
                                        ! to be changed (check it one day)
!-----------------------------------[Locals]-----------------------------------!
  type(Comm_Type), pointer :: Comm
  character(SL)            :: vn
  integer                  :: vs, disp_loop, cnt_loop, nf
  integer                  :: s, c1, c2, c1g, c2g, cnt, i_sid, sg
  real                     :: buffer
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
      if(First_Proc()) print *, '# Reading variable: ', trim(vn)

      i_sid = 1                 ! local side counter
      do sg = 1, Comm % nf_tot  ! global side loop

        ! The global value is read by all processors
        call Comm % Read_Real(fh, buffer, disp_loop)

        ! If you didn't go beyond local face map ...
        ! ... check if this global side is present in this processor
        if(i_sid .le. Grid % n_faces) then
          if(Comm % face_map_dup_glo(i_sid) .eq. sg) then
            s = Comm % face_map_dup_loc(i_sid)  ! get the local face number
            c1 = Grid % faces_c(1,s)
            c2 = Grid % faces_c(2,s)
            c1g = Grid % Comm % cell_glo(c1)
            c2g = Grid % Comm % cell_glo(c2)
            array(s) = buffer
            if(c2g > 0 .and. c2g < c1g) array(s) = -buffer
            i_sid = i_sid + 1
          end if
        end if

      end do

      disp = disp_loop

      return

    ! If variable not found, advance the offset only
    else
      disp_loop = disp_loop + vs
    end if

    ! Check if variable is in the file
    if(cnt_loop >= vc) goto 1

  end do

  return

1 if(First_Proc()) print *, '# Variable: ', trim(var_name), ' not found! ',  &
                             'Setting the values to 0.0!'
  array(:) = 0.0

  end subroutine
