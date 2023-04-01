!==============================================================================!
  subroutine Iter_Fill_Scalar_At(Info, r, c, name_var, res, n_iter)
!------------------------------------------------------------------------------!
!   Inserts infromation at specified position in the information box.          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Info_Type), intent(out) :: Info
  integer,          intent(in)  :: r         ! row
  integer,          intent(in)  :: c         ! column
  character(len=*), intent(in)  :: name_var
  real,             intent(in)  :: res
  integer,          intent(in)  :: n_iter    ! number of iterations
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
