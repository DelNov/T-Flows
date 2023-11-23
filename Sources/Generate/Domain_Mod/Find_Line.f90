!==============================================================================!
  subroutine Find_Line(Dom, n1, n2, res)
!------------------------------------------------------------------------------!
!>  Searches for the smallest block in domain (Dom) where the line defined by
!>  nodes n1 and n2 exists, and outputs its resolution (res).
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Domain_Type)   :: Dom     !! computational domain
  integer, intent(in)  :: n1, n2  !! node
  integer, intent(out) :: res     !! output resolution
!-----------------------------------[Locals]-----------------------------------!
  integer :: b, l1, l2
!==============================================================================!

  do b = 1, size(Dom % blocks)
    do l1 = 1, 8
      do l2 = 1, 8
        if( (Dom % blocks(b) % corners(l1) .eq. n1) .and. &
            (Dom % blocks(b) % corners(l2) .eq. n2) ) then
          if( iabs(l2-l1) .eq. 1 ) res = Dom % blocks(b) % resolutions(1)
          if( iabs(l2-l1) .eq. 2 ) res = Dom % blocks(b) % resolutions(2)
          if( iabs(l2-l1) .eq. 4 ) res = Dom % blocks(b) % resolutions(3)
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
