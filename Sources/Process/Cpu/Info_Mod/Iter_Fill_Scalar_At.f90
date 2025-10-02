!==============================================================================!
  subroutine Iter_Fill_Scalar_At(Info, r, c, name_var, res, n_iter)
!------------------------------------------------------------------------------!
!>  Iter_Fill_Scalar_At is designed to insert detailed information about
!>  user-defined scalar variables at a specified position within the iteration
!>  information box. It updates the contents of the box with the name of the
!>  scalar variable, its residual value, and the number of iterations it has
!>  undergone.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Info_Type)             :: Info       !! parent, singleton object Info
  integer,          intent(in) :: r          !! row in the info box
  integer,          intent(in) :: c          !! column in the info box
  character(len=*),  intent(in) :: name_var  !! name of the variable to print
  real,              intent(in) :: res       !! residual after linear solver
  integer, optional, intent(in) :: n_iter    !! iterations performed
                                             !! in the linear solver
!==============================================================================!

  if(First_Proc()) then

    ! Update the number of lines you'll have to print in the end
    Info % iter % n_user_lines = max(Info % iter % n_user_lines, r)

    ! User variables
    write(Info % iter % lines_user(r)((c-1)*L_BOX+ 3 :  &
                                      (c-1)*L_BOX+ 6),  '(a4)')  name_var
    write(Info % iter % lines_user(r)((c-1)*L_BOX+ 7 :  &
                                      (c-1)*L_BOX+ 7),  '(a1)')  ':'
    write(Info % iter % lines_user(r)((c-1)*L_BOX+ 8 :  &
                                      (c-1)*L_BOX+10),  '(i3)')  n_iter
    ! Residual
    write(Info % iter % lines_user(r)((c-1)*L_BOX+12 :  &
                                      (c-1)*L_BOX+20),  '(1pe9.3)') res

  end if

  end subroutine
