!==============================================================================!
  subroutine Info_Mod_Iter_Fill_User_At(r, c, name_var, n_iter, res)
!------------------------------------------------------------------------------!
!   Inserts infromation at specified position in the information box.          !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Comm_Mod, only: this_proc    
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer          :: r         ! row
  integer          :: c         ! column
  character(len=*) :: name_var
  integer          :: n_iter    ! number of iterations
  real             :: res       
!==============================================================================!

  if (this_proc < 2) then

    ! Update the number of lines you'll have to print in the end
    iter_info % n_user_lines = max(iter_info % n_user_lines, r)

    ! User variables
    write(iter_info % lines_user(r)((c-1)*L_BOX+ 3 :  &
                                    (c-1)*L_BOX+ 6),  '(a4)')  name_var
    write(iter_info % lines_user(r)((c-1)*L_BOX+ 7 :  &
                                    (c-1)*L_BOX+ 7),  '(a1)')  ':'
    write(iter_info % lines_user(r)((c-1)*L_BOX+ 8 :  &
                                    (c-1)*L_BOX+10),  '(i3)')  n_iter
    ! Residual 
    write(iter_info % lines_user(r)((c-1)*L_BOX+12 :  &
                                    (c-1)*L_BOX+20),  '(1pe9.3)') res

  end if

  end subroutine
