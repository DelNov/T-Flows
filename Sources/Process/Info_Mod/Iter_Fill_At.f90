!==============================================================================!
  subroutine Info_Mod_Iter_Fill_At(r, c, name_var, n_iter, res)
!------------------------------------------------------------------------------!
!   Inserts infromation at specified position in the information box.          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer          :: r         ! row
  integer          :: c         ! column
  character(len=*) :: name_var
  integer          :: n_iter    ! number of iterations
  real             :: res
!==============================================================================!

  if(this_proc < 2) then

    ! Normal variables
    if(n_iter > -1) then  ! if n_iter .eq. -1, iterations won't be written
                          ! (that's useful for mass residual)
      write(iter_info % lines(r)((c-1)*L_BOX+ 3 :  &
                                 (c-1)*L_BOX+ 6),  '(a4)')  name_var
      write(iter_info % lines(r)((c-1)*L_BOX+ 7 :  &
                                 (c-1)*L_BOX+ 7),  '(a1)')  ':'
      write(iter_info % lines(r)((c-1)*L_BOX+ 8 :  &
                                 (c-1)*L_BOX+10),  '(i3)')  n_iter
    endif

    ! Residual 
    write(iter_info % lines(r)((c-1)*L_BOX+12 :  &
                               (c-1)*L_BOX+20),  '(1pe9.3)') res

  end if

  end subroutine
