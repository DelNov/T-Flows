!==============================================================================!
  subroutine Info_Mod_Iter_Start()
!------------------------------------------------------------------------------!
!   Essentially creates a box in which iteration residuls will be written.     !
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  integer :: i, l
!==============================================================================!

  if(First_Proc()) then

    ! Create frame
    do i = 1, L_LINE, L_LINE-1
      iter_info % line_lead      (i:i) = '#'
      iter_info % line_iter      (i:i) = '#'
      iter_info % line_sep       (i:i) = '#'
    end do

    ! For normal lines
    do l = 1, 4
      iter_info % lines(l)(     1:L_LINE) = ' '
      iter_info % lines(l)(     1:     1) = '#'
      iter_info % lines(l)(L_LINE:L_LINE) = '#'
    end do

    ! For user lines
    do l = 1, MAX_USER_LINES
      iter_info % lines_user(l)(     1:L_LINE) = ' '
      iter_info % lines_user(l)(     1:     1) = '#'
      iter_info % lines_user(l)(L_LINE:L_LINE) = '#'
    end do

    ! Lead and separating lines
    do i = 2, L_LINE-1
      iter_info % line_lead(i:i) = '='
      iter_info % line_sep (i:i) = '-'
    end do

    ! Create separators (character must be length of L_BOX)
    do i = 2, L_LINE-L_BOX, L_BOX

      ! For normal lines
      do l = 1, 4
        write(iter_info % lines(l)      (i:i+L_BOX-1), '(a21)') '|'
      end do

      ! For user lines
      do l = 1, MAX_USER_LINES
        write(iter_info % lines_user(l) (i:i+L_BOX-1), '(a21)') '|'
      end do
      write(iter_info % line_sep(i+L_BOX-1 :  &
                                 i+L_BOX-1),   '(a1)') '+'
    end do

  end if

  end subroutine
