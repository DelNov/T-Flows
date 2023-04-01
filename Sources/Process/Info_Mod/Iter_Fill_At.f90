!==============================================================================!
  subroutine Iter_Fill_At(Info, r, c, name_var, res, n_iter)
!------------------------------------------------------------------------------!
!   Inserts infromation at specified position in the information box.          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Info_Type),  intent(out) :: Info
  integer,           intent(in)  :: r         ! row
  integer,           intent(in)  :: c         ! column
  character(len=*),  intent(in)  :: name_var
  real,              intent(in)  :: res
  integer, optional, intent(in)  :: n_iter    ! number of iterations
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
