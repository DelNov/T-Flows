!==============================================================================!
  subroutine Domain_Mod_Find_Line(dom, n1, n2, res)
!------------------------------------------------------------------------------!
!   Searches for a smallest block where the line defined by n1-n2 is.          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Domain_Type)    :: dom
  integer, intent(in)  :: n1, n2
  integer, intent(out) :: res
!-----------------------------------[Locals]-----------------------------------!
  integer :: b, l1, l2
!==============================================================================!

  do b = 1, size(dom % blocks)
    do l1 = 1, 8
      do l2 = 1, 8
        if( (dom % blocks(b) % corners(l1) .eq. n1) .and. &
            (dom % blocks(b) % corners(l2) .eq. n2) ) then
          if( iabs(l2-l1) .eq. 1 ) res = dom % blocks(b) % resolutions(1)
          if( iabs(l2-l1) .eq. 2 ) res = dom % blocks(b) % resolutions(2)
          if( iabs(l2-l1) .eq. 4 ) res = dom % blocks(b) % resolutions(3)
          goto 1
        end if
      end do
    end do
  end do

  print *, '# ERROR in Generator'
  print *, '# You tried to define the line', n1, n2, ' but it'
  print *, '# doesn''t exists in the block specifications.'
  print *, '# Exiting !'
  stop

1 return

  end subroutine
