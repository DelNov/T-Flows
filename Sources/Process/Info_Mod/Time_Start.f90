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
    time_info % line_lead (2*L_BOX-2 : 2*L_BOX-2) = '#'
    time_info % lines(1)  (2*L_BOX-2 : 2*L_BOX-2) = '#'
    time_info % lines(2)  (2*L_BOX-2 : 2*L_BOX-2) = '#'
    time_info % lines(3)  (2*L_BOX-2 : 2*L_BOX-2) = '#'
    time_info % lines(4)  (2*L_BOX-2 : 2*L_BOX-2) = '#'
    time_info % lines(5)  (2*L_BOX-2 : 2*L_BOX-2) = '#'
    time_info % lines(6)  (2*L_BOX-2 : 2*L_BOX-2) = '#'
    time_info % line_trail(2*L_BOX-2 : 2*L_BOX-2) = '#'

    time_info % line_lead (4*L_BOX+4 : 4*L_BOX+4) = '#'
    time_info % lines(1)  (4*L_BOX+4 : 4*L_BOX+4) = '#'
    time_info % lines(2)  (4*L_BOX+4 : 4*L_BOX+4) = '#'
    time_info % lines(3)  (4*L_BOX+4 : 4*L_BOX+4) = '#'
    time_info % lines(4)  (4*L_BOX+4 : 4*L_BOX+4) = '#'
    time_info % lines(5)  (4*L_BOX+4 : 4*L_BOX+4) = '#'
    time_info % lines(6)  (4*L_BOX+4 : 4*L_BOX+4) = '#'
    time_info % line_trail(4*L_BOX+4 : 4*L_BOX+4) = '#'

    do i = 2*L_BOX-1, 4*L_BOX+3
      time_info % line_lead (i:i) = '='
      time_info % line_trail(i:i) = '-'
    end do

  end if

  end subroutine
