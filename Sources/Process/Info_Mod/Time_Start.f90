!==============================================================================!
  subroutine Time_Start(Info)
!------------------------------------------------------------------------------!
!   Essentially creates a box in which iteration residuls will be written.     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Info_Type) :: Info
!-----------------------------------[Locals]-----------------------------------!
  integer :: i
!==============================================================================!

  if(First_Proc()) then

    ! Create a frame
    Info % time % line_lead (2*L_BOX-2 : 2*L_BOX-2) = '#'
    Info % time % line(1)   (2*L_BOX-2 : 2*L_BOX-2) = '#'
    Info % time % line(2)   (2*L_BOX-2 : 2*L_BOX-2) = '#'
    Info % time % line(3)   (2*L_BOX-2 : 2*L_BOX-2) = '#'
    Info % time % line(4)   (2*L_BOX-2 : 2*L_BOX-2) = '#'
    Info % time % line(5)   (2*L_BOX-2 : 2*L_BOX-2) = '#'
    Info % time % line(6)   (2*L_BOX-2 : 2*L_BOX-2) = '#'
    Info % time % line_trail(2*L_BOX-2 : 2*L_BOX-2) = '#'

    Info % time % line_lead (4*L_BOX+4 : 4*L_BOX+4) = '#'
    Info % time % line(1)   (4*L_BOX+4 : 4*L_BOX+4) = '#'
    Info % time % line(2)   (4*L_BOX+4 : 4*L_BOX+4) = '#'
    Info % time % line(3)   (4*L_BOX+4 : 4*L_BOX+4) = '#'
    Info % time % line(4)   (4*L_BOX+4 : 4*L_BOX+4) = '#'
    Info % time % line(5)   (4*L_BOX+4 : 4*L_BOX+4) = '#'
    Info % time % line(6)   (4*L_BOX+4 : 4*L_BOX+4) = '#'
    Info % time % line_trail(4*L_BOX+4 : 4*L_BOX+4) = '#'

    do i = 2*L_BOX-1, 4*L_BOX+3
      Info % time % line_lead (i:i) = '='
      Info % time % line_trail(i:i) = '-'
    end do

  end if

  end subroutine
