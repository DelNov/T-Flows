!==============================================================================!
  subroutine Backup_Mod_Write_Face_Real(Grid, fh, disp, vc, var_name, array, &
                                        correct_sign)
!------------------------------------------------------------------------------!
!   Writes a vector variable with face values backup file.                     !
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
  integer                  :: vs  ! variable size
  integer                  :: s, c1, c2, cg1, cg2, sg, i_sid, error
  real                     :: buffer
!==============================================================================!

  ! Take alias
  Comm => Grid % Comm

  if(this_proc < 2) print *, '# Writing variable: ', trim(var_name)

  ! Increase variable count
  vc = vc + 1

  ! Vector with boundaries
  vn = var_name
  call Comm % Write_Text(fh, vn, disp)

  vs = (Comm % nf_tot) * SIZE_REAL
  call Comm % Write_Int (fh, vs, disp)

  i_sid = 1
  do sg = 1, Comm % nf_tot
    buffer = 0.0

    ! If you didn't go beyond local face map ...
    ! ... check if this global side is present in this processor
    if(i_sid .le. Comm % nf_sub) then
      if(Comm % face_map_uni_glo(i_sid) .eq. sg) then
        s = Comm % face_map_uni_loc(i_sid)  ! get the local face number
        buffer = array(s)
        i_sid = i_sid + 1
      end if
    end if
    call Comm_Mod_Global_Sum_Real(buffer)     ! distribute the value
    call Comm % Write_Real(fh, buffer, disp)  ! write the same value from all
  end do

  end subroutine
