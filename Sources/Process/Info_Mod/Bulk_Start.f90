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
      Info % bulk % line_lead(i:i) = '-'
      Info % bulk % line_foll(i:i) = ' '
    end do
    Info % bulk % line_lead(     1:     1) = '#'
    Info % bulk % line_lead(L_LINE:L_LINE) = '#'

    do i = L_BOX + 1,  5*L_BOX + 1
      Info % bulk % line_lead (i:i) = '='
      Info % bulk % line_foll (i:i) = '='
      Info % bulk % line_sep  (i:i) = '-'
      Info % bulk % line_trail(i:i) = '-'
    end do

    Info % bulk % line_lead(  L_BOX+1:  L_BOX+1) = '+'
    Info % bulk % line_foll(  L_BOX+1:  L_BOX+1) = '+'
    Info % bulk % line_lead(5*L_BOX+1:5*L_BOX+1) = '+'
    Info % bulk % line_foll(5*L_BOX+1:5*L_BOX+1) = '+'

    do i = 22, 106, 106-22
      Info % bulk % line(1)   (i:i) = '#'
      Info % bulk % line_sep  (i:i) = '#'
      Info % bulk % line(2)   (i:i) = '#'
      Info % bulk % line(3)   (i:i) = '#'
      Info % bulk % line_trail(i:i) = '#'
    end do

    ! Create separators
    Info % bulk % line_lead(L_LINE/2+1:L_LINE/2+1) = '+'
    Info % bulk % line_foll(L_LINE/2+1:L_LINE/2+1) = '+'
    Info % bulk % line(1)  (L_LINE/2+1:L_LINE/2+1) = '|'
    Info % bulk % line_sep (L_LINE/2+1:L_LINE/2+1) = '+'

    Info % bulk % line_sep  (50:50) = '+'
    Info % bulk % line(2)   (50:50) = '|'
    Info % bulk % line(3)   (50:50) = '|'
    Info % bulk % line_trail(50:50) = '+'

    Info % bulk % line_sep  (78:78) = '+'
    Info % bulk % line(2)   (78:78) = '|'
    Info % bulk % line(3)   (78:78) = '|'
    Info % bulk % line_trail(78:78) = '+'

  end if

  end subroutine
