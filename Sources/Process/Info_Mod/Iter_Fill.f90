!==============================================================================!
  subroutine Info_Mod_Iter_Fill(n)
!------------------------------------------------------------------------------!
!   Fills the info box with information to be written on the screen.           !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Comm_Mod, only: this_proc    
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: n  ! inner iteration
!-----------------------------------[Locals]-----------------------------------!
  integer :: c
!==============================================================================!

  if (this_proc < 2) then

    ! Write basic info
    c = 1  ! a column
    write(iter_info % lines(1)((c-1)*L_BOX+ 5 :  &
                               (c-1)*L_BOX+14),  '(a10)') 'Iteration:'
    write(iter_info % lines(1)((c-1)*L_BOX+15 :  &
                               (c-1)*L_BOX+17),   '(i3)') n

    c = 2  ! a column
    write(iter_info % lines(1)((c-1)*L_BOX+3 :  &
                               (c-1)*L_BOX+7),    '(a5)') 'Mass:'

  end if

  end subroutine
