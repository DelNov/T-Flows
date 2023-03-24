!==============================================================================!
  subroutine Info_Mod_Iter_Fill(n)
!------------------------------------------------------------------------------!
!   Fills the info box with information to be written on the screen.           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: n  ! inner iteration
!-----------------------------------[Locals]-----------------------------------!
  integer :: c
!==============================================================================!

  if(First_Proc()) then

    ! Write basic info
    c = 1  ! a column
    write(iter_info % line_iter((c-1)*L_BOX+58 :  &
                                (c-1)*L_BOX+67),  '(a10)') 'Iteration:'
    write(iter_info % line_iter((c-1)*L_BOX+68 :  &
                                (c-1)*L_BOX+70),   '(i3)') n

    c = 5  ! a column
    write(iter_info % lines(1)((c-1)*L_BOX+3 :  &
                               (c-1)*L_BOX+7),    '(a5)') 'MASS:'

  end if

  end subroutine
