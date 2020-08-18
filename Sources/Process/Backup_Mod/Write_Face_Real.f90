!==============================================================================!
  subroutine Backup_Mod_Write_Face_Real(grid, fh, disp, vc, var_name, array, &
                                        correct_sign)
!------------------------------------------------------------------------------!
!   Writes a vector variable with face values backup file.                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type), target :: grid
  integer                 :: fh, disp, vc
  character(len=*)        :: var_name
  real                    :: array(grid % n_faces)
  logical, optional       :: correct_sign  ! for face fluxes, signs might have
                                           ! to be changed (check it one day)
!-----------------------------------[Locals]-----------------------------------!
  type(Comm_Type), pointer :: comm
  character(len=80)        :: vn
  integer                  :: vs  ! variable size
  integer                  :: s, c1, c2, cg1, cg2
!==============================================================================!

  ! Take alias
  comm => grid % comm

  if(this_proc < 2) print *, '# Writing variable: ', trim(var_name)

  ! Correct the signs of fluxes  (Remember, they are defined
  ! defined to be positive from cg1 to cg2; and cg2 > cg1)
  if(present(correct_sign)) then
    if(correct_sign) then
      do s = 1, grid % n_faces
        c1  = grid % faces_c(1,s)
        c2  = grid % faces_c(2,s)
        cg1 = grid % comm % cell_glo(c1)
        cg2 = grid % comm % cell_glo(c2)
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
  call Comm_Mod_Write_Text(fh, vn, disp)

  vs = (comm % nf_t) * SIZE_REAL
  call Comm_Mod_Write_Int (fh, vs, disp)

  call Comm_Mod_Write_Face_Real(comm, fh, array(1:comm % nf_s), disp)

  ! Correct the signs of fluxes  (Remember, they are defined
  ! defined to be positive from cg1 to cg2; and cg2 > cg1)
  if(present(correct_sign)) then
    if(correct_sign) then
      do s = 1, grid % n_faces
        c1  = grid % faces_c(1,s)
        c2  = grid % faces_c(2,s)
        cg1 = grid % comm % cell_glo(c1)
        cg2 = grid % comm % cell_glo(c2)
        if(cg2 > 0 .and. cg2 < cg1) then
          array(s) = -array(s)
        end if
      end do
    end if
  end if

  end subroutine
