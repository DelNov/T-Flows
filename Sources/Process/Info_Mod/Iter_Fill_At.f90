!==============================================================================!
  subroutine Iter_Fill_At(Info, r, c, name_var, res, n_iter)
!------------------------------------------------------------------------------!
!>  The Iter_Fill_At subroutine populates a specific row and column of the
!>  iteration information box with details such as the name of the variable,
!>  its residual, and the number of iterations (if provided). This function is
!>  vital for displaying specific simulation data within the formatted
!>  iteration box. It ensures that iteration-related metrics are organized and
!>  presented in a clear, concise manner.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Info_Type)              :: Info      !! parent, singleton object Info
  integer,           intent(in) :: r         !! row in the info box
  integer,           intent(in) :: c         !! column in the info box
  character(len=*),  intent(in) :: name_var  !! name of the variable to print
  real,              intent(in) :: res       !! residual after linear solver
  integer, optional, intent(in) :: n_iter    !! iterations performed
                                             !! in the linear solver
!==============================================================================!

  if(First_Proc()) then

    ! Normal variables
    if(present(n_iter)) then  ! if n_iter not present, iterations won't be
                              ! written (that's useful for mass residual)
      write(Info % iter % line(r)((c-1)*L_BOX+ 3 :  &
                                  (c-1)*L_BOX+ 6),  '(a4)')  name_var
      write(Info % iter % line(r)((c-1)*L_BOX+ 7 :  &
                                  (c-1)*L_BOX+ 7),  '(a1)')  ':'
      write(Info % iter % line(r)((c-1)*L_BOX+ 8 :  &
                                  (c-1)*L_BOX+10),  '(i3)')  n_iter
    endif

    ! Residual
    write(Info % iter % line(r)((c-1)*L_BOX+12 :  &
                                (c-1)*L_BOX+20),  '(1pe9.3)') res

  end if

  end subroutine
