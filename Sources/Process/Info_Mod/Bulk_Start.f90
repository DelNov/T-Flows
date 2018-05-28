!==============================================================================!
  subroutine Info_Mod_Bulk_Start()
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
    do i = 1, L_LINE, L_LINE-1
      bulk_info % line_lead (i:i) = '#'
      bulk_info % lines(1)  (i:i) = '#'
      bulk_info % line_sep  (i:i) = '#'
      bulk_info % lines(2)  (i:i) = '#'
      bulk_info % lines(3)  (i:i) = '#'
      bulk_info % line_trail(i:i) = '#'
    end do

    do i = 2, L_LINE-1
      bulk_info % line_lead (i:i) = '='
      bulk_info % line_sep  (i:i) = '-'
      bulk_info % line_trail(i:i) = '-'
    end do

    ! Create separators
    bulk_info % line_lead(L_LINE/2+1:L_LINE/2+1) = '+'
    bulk_info % lines(1) (L_LINE/2+1:L_LINE/2+1) = '|'
    bulk_info % line_sep (L_LINE/2+1:L_LINE/2+1) = '+'

    bulk_info % line_sep  (29:29) = '+'
    bulk_info % lines(2)  (29:29) = '|'
    bulk_info % lines(3)  (29:29) = '|'
    bulk_info % line_trail(29:29) = '+'

    bulk_info % line_sep  (57:57) = '+'
    bulk_info % lines(2)  (57:57) = '|'
    bulk_info % lines(3)  (57:57) = '|'
    bulk_info % line_trail(57:57) = '+'

  end if

  end subroutine
