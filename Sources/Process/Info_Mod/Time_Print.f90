!==============================================================================!
  subroutine Info_Mod_Time_Print()
!------------------------------------------------------------------------------!
!   Prints information about inner iteration on the screen.                    !
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  integer :: i
!==============================================================================!

  if(this_proc < 2) then

    print '(a129)', ' '
    print '(a129)', ' '
    print '(a129)', time_info % line_lead

    ! Print only lines which have colon in the first column :-)
    do i=1,6
      print '(a129)', time_info % lines(i)
    end do

    print '(a129)', time_info % line_trail
    print '(a129)', ' '

  end if

  end subroutine
