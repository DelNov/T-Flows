!==============================================================================!
  subroutine Iter_Start(Info)
!------------------------------------------------------------------------------!
!   Essentially creates a box in which iteration residuls will be written.     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Info_Type), intent(out) :: Info
!-----------------------------------[Locals]-----------------------------------!
  integer :: i, l
!==============================================================================!

  if(First_Proc()) then

    ! Create frame
    do i = 1, L_LINE, L_LINE-1
      Info % iter % line_lead      (i:i) = '#'
      Info % iter % line_iter      (i:i) = '#'
      Info % iter % line_sep       (i:i) = '#'
    end do

    ! For normal lines
    do l = 1, 4
      Info % iter % line(l)(     1:L_LINE) = ' '
      Info % iter % line(l)(     1:     1) = '#'
      Info % iter % line(l)(L_LINE:L_LINE) = '#'
    end do

    ! For user lines
    do l = 1, MAX_USER_LINES
      Info % iter % lines_user(l)(     1:L_LINE) = ' '
      Info % iter % lines_user(l)(     1:     1) = '#'
      Info % iter % lines_user(l)(L_LINE:L_LINE) = '#'
    end do

    ! Lead and separating lines
    do i = 2, L_LINE-1
      Info % iter % line_lead(i:i) = '='
      Info % iter % line_sep (i:i) = '-'
    end do

    ! Create separators (character must be length of L_BOX)
    do i = 2, L_LINE-L_BOX, L_BOX

      ! For normal lines
      do l = 1, 4
        write(Info % iter % line(l)       (i:i+L_BOX-1), '(a21)') '|'
      end do

      ! For user lines
      do l = 1, MAX_USER_LINES
        write(Info % iter % lines_user(l) (i:i+L_BOX-1), '(a21)') '|'
      end do
      write(Info % iter % line_sep(i+L_BOX-1 :  &
                                   i+L_BOX-1),   '(a1)') '+'
    end do

  end if

  end subroutine
