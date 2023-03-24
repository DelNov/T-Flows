!==============================================================================!
  subroutine Info_Mod_Bulk_Start()
!------------------------------------------------------------------------------!
!   Essentially creates a box in which iteration residuls will be written.     !
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  integer :: i
!==============================================================================!

  if(First_Proc()) then

    ! Create a frame

    do i = 1, L_LINE
      bulk_info % line_lead(i:i) = '-'
      bulk_info % line_foll(i:i) = ' '
    end do
    bulk_info % line_lead(     1:     1) = '#'
    bulk_info % line_lead(L_LINE:L_LINE) = '#'

    do i = L_BOX + 1,  5*L_BOX + 1
      bulk_info % line_lead (i:i) = '='
      bulk_info % line_foll (i:i) = '='
      bulk_info % line_sep  (i:i) = '-'
      bulk_info % line_trail(i:i) = '-'
    end do

    bulk_info % line_lead(  L_BOX+1:  L_BOX+1) = '+'
    bulk_info % line_foll(  L_BOX+1:  L_BOX+1) = '+'
    bulk_info % line_lead(5*L_BOX+1:5*L_BOX+1) = '+'
    bulk_info % line_foll(5*L_BOX+1:5*L_BOX+1) = '+'

    do i = 22, 106, 106-22
      bulk_info % lines(1)  (i:i) = '#'
      bulk_info % line_sep  (i:i) = '#'
      bulk_info % lines(2)  (i:i) = '#'
      bulk_info % lines(3)  (i:i) = '#'
      bulk_info % line_trail(i:i) = '#'
    end do

    ! Create separators
    bulk_info % line_lead(L_LINE/2+1:L_LINE/2+1) = '+'
    bulk_info % line_foll(L_LINE/2+1:L_LINE/2+1) = '+'
    bulk_info % lines(1) (L_LINE/2+1:L_LINE/2+1) = '|'
    bulk_info % line_sep (L_LINE/2+1:L_LINE/2+1) = '+'

    bulk_info % line_sep  (50:50) = '+'
    bulk_info % lines(2)  (50:50) = '|'
    bulk_info % lines(3)  (50:50) = '|'
    bulk_info % line_trail(50:50) = '+'

    bulk_info % line_sep  (78:78) = '+'
    bulk_info % lines(2)  (78:78) = '|'
    bulk_info % lines(3)  (78:78) = '|'
    bulk_info % line_trail(78:78) = '+'

  end if

  end subroutine
