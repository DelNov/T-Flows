!==============================================================================!
  subroutine Load_Shadow_Faces(grid)
!------------------------------------------------------------------------------!
! Reads:  name.shadow                                                          !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Name_Mod, only: problem_name
  use gen_mod 
  use Div_Mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer           :: c1, c2, n, s, dum_i
  character(len=80) :: dum_s, name_in
!==============================================================================!

  !----------------------!    
  !   Read shadow file   !
  !----------------------!    
  name_in = problem_name
  call Name_File(0, name_in, '.shadow')
  open(9, file=name_in)
  print *, '# Reading the file: ', name_in

  do s = grid % n_faces+1,grid % n_faces + grid % n_sh
    read(9,*) grid % faces_n_nodes(s)
    if(grid % faces_n_nodes(s) .eq. 3) then
      read(9,*) grid % faces_n(1,s),  &
                grid % faces_n(2,s),  &
                grid % faces_n(3,s),  &
                grid % faces_c(1,s), grid % faces_c(2,s)
    else if(grid % faces_n_nodes(s) .eq. 4) then
      read(9,*) grid % faces_n(1,s),  &
                grid % faces_n(2,s),  &
                grid % faces_n(3,s),  &
                grid % faces_n(4,s),  &
                grid % faces_c(1,s), grid % faces_c(2,s)
    end if
  end do

  close(9)

  end subroutine
