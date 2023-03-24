!==============================================================================!
  subroutine Info_Mod_Time_Print()
!------------------------------------------------------------------------------!
!   Prints information about inner iteration on the screen.                    !
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  integer :: i
!==============================================================================!

  if(First_Proc()) then

    print '(a)',   ''
    print '(a)',   ''
    print '(a90)', trim(time_info % line_lead)

    ! Print only lines which have colon in the first column :-)
    do i=1,6
      print '(a90)', trim(time_info % lines(i))
    end do

    print '(a90)', trim(time_info % line_trail)
    print '(a)',   ''

  end if

  end subroutine
