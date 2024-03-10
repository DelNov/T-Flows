!==============================================================================!
  subroutine Time_Print(Info)
!------------------------------------------------------------------------------!
!>  Responsible for printing the time information box, previously filled with
!>  current simulation time data, on the screen. It ensures that the
!>  time-related data is visually accessible and readable during a simulation.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Info_Type) :: Info  !! parent, singleton object Info
!-----------------------------------[Locals]-----------------------------------!
  integer :: i
!==============================================================================!

  if(First_Proc()) then

    print '(a)',   ''
    print '(a)',   ''
    print '(a90)', trim(Info % time % line_lead)

    ! Print only lines which have colon in the first column :-)
    do i = 1, 6
      print '(a90)', trim(Info % time % line(i))
    end do

    print '(a90)', trim(Info % time % line_trail)
    print '(a)',   ''

  end if

  end subroutine
