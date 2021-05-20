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
  integer                  :: s, c1, c2, cg1, cg2
!==============================================================================!

  ! Take alias
  Comm => Grid % Comm

  if(this_proc < 2) print *, '# Writing variable: ', trim(var_name)

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

  ! Increase variable count
  vc = vc + 1

  ! Vector with boundaries
  vn = var_name
  call Comm % Write_Text(fh, vn, disp)

  vs = (Comm % nf_tot) * SIZE_REAL
  call Comm % Write_Int (fh, vs, disp)

  call Comm % Write_Face_Real(fh, array(1:Comm % nf_sub), disp)

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

  end subroutine
