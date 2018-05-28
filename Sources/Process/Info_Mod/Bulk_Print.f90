!==============================================================================!
  subroutine Info_Mod_Bulk_Print()
!------------------------------------------------------------------------------!
!   Prints information about inner iteration on the screen.                    !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Comm_Mod, only: this_proc    
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  if (this_proc < 2) then

    ! String is L_LINE+2 long
    print '(a87)', bulk_info % line_lead  
    print '(a87)', bulk_info % lines(1)
    print '(a87)', bulk_info % line_sep
    print '(a87)', bulk_info % lines(2)
    print '(a87)', bulk_info % lines(3)
    print '(a87)', bulk_info % line_trail  
    print '(a87)', ' '

  end if
                 
  end subroutine
