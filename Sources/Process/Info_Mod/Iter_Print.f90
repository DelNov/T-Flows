!==============================================================================!
  subroutine Iter_Print(Info, d)
!------------------------------------------------------------------------------!
!   Prints information about inner iteration on the screen.                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Info_Type)    :: Info
  integer, intent(in) :: d  ! domain
!-----------------------------------[Locals]-----------------------------------!
  integer               :: i
  character(len=L_LINE) :: tmp
!==============================================================================!

  if(First_Proc()) then

    if(d .eq. 1) then
      print '(a129)', Info % iter % line_lead
      print '(a129)', Info % iter % line_iter
      print '(a129)', Info % iter % line_sep
    end if

    ! Print only lines which have colon in the first column :-)
    print '(a129)', Info % iter % line(1)

    ! First print normal lines
    do i = 2, 4
      tmp = Info % iter % line(i)
      if( tmp(7:7) .eq. ':') print '(a129)', Info % iter % line(i)
      Info % iter % line(i)(7:7) = ' '        ! remove column for other domains
    end do

    ! Then print user lines
    do i = 1, Info % iter % n_user_lines
      tmp = Info % iter % lines_user(i)
      if( tmp(7:7) .eq. ':') print '(a129)', Info % iter % lines_user(i)
      Info % iter % lines_user(i)(7:7) = ' '  ! remove column for other domains
    end do

    call Info % Iter_Start()

  end if

  end subroutine
