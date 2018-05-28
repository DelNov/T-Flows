!==============================================================================!
  subroutine Info_Mod_Time_Start()
!------------------------------------------------------------------------------!
!   Essentially creates a box in which iteration residuls will be written.     !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Comm_Mod, only: this_proc
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  integer :: i
!==============================================================================!

  if (this_proc < 2) then

    ! Create a frame
    time_info % line_lead (L_BOX+1:L_BOX+1) = '#'
    time_info % lines(1)  (L_BOX+1:L_BOX+1) = '#'
    time_info % lines(2)  (L_BOX+1:L_BOX+1) = '#'
    time_info % lines(3)  (L_BOX+1:L_BOX+1) = '#'
    time_info % lines(4)  (L_BOX+1:L_BOX+1) = '#'
    time_info % lines(5)  (L_BOX+1:L_BOX+1) = '#'
    time_info % lines(6)  (L_BOX+1:L_BOX+1) = '#'
    time_info % line_trail(L_BOX+1:L_BOX+1) = '#'

    time_info % line_lead (L_LINE-L_BOX:L_LINE-L_BOX) = '#'
    time_info % lines(1)  (L_LINE-L_BOX:L_LINE-L_BOX) = '#'
    time_info % lines(2)  (L_LINE-L_BOX:L_LINE-L_BOX) = '#'
    time_info % lines(3)  (L_LINE-L_BOX:L_LINE-L_BOX) = '#'
    time_info % lines(4)  (L_LINE-L_BOX:L_LINE-L_BOX) = '#'
    time_info % lines(5)  (L_LINE-L_BOX:L_LINE-L_BOX) = '#'
    time_info % lines(6)  (L_LINE-L_BOX:L_LINE-L_BOX) = '#'
    time_info % line_trail(L_LINE-L_BOX:L_LINE-L_BOX) = '#'

    do i = L_BOX+2, L_LINE-L_BOX-1
      time_info % line_lead (i:i) = '='
      time_info % line_trail(i:i) = '-'
    end do

  end if

  end subroutine
