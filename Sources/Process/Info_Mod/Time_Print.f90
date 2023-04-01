!==============================================================================!
  subroutine Time_Print(Info)
!------------------------------------------------------------------------------!
!   Prints information about inner iteration on the screen.                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Info_Type), intent(out) :: Info
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
