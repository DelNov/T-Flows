!==============================================================================!
  subroutine Info_Mod_Time_Print()
!------------------------------------------------------------------------------!
!   Prints information about inner iteration on the screen.                    !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Comm_Mod, only: this_proc    
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  integer               :: i
!==============================================================================!

  if (this_proc < 2) then

    print '(a87)', ' '
    print '(a87)', time_info % line_lead  

    ! Print only lines which have colon in the first column :-)
    do i=1,6
      print '(a87)', time_info % lines(i)
    end do

    print '(a87)', time_info % line_trail  
    print '(a87)', ' '

  end if
                 
  end subroutine
