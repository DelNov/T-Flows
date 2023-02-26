!==============================================================================!
  subroutine Find_Surface(Dom, n1, n2, n3, n4, block, face)
!------------------------------------------------------------------------------!
!   Searches for a block where the surface defined by n1, n2, n3, n4 is.       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Domain_Type)   :: Dom
  integer, intent(in)  :: n1, n2, n3, n4
  integer, intent(out) :: block, face
!-----------------------------------[Locals]-----------------------------------!
  integer :: b, fc, p1, p2, p3, p4
!==============================================================================!

  do b = 1, size(Dom % blocks)
    do fc = 1, 6
      p1 = Dom % blocks(b) % faces(fc, 1)
      p2 = Dom % blocks(b) % faces(fc, 2)
      p3 = Dom % blocks(b) % faces(fc, 3)
      p4 = Dom % blocks(b) % faces(fc, 4)
      if( ((p1 .eq. n1).and.(p3 .eq. n3)) .or.  &
          ((p1 .eq. n4).and.(p3 .eq. n2)) .or.  &
          ((p1 .eq. n3).and.(p3 .eq. n1)) .or.  &
          ((p1 .eq. n2).and.(p3 .eq. n4)) ) goto 1
    end do
  end do 

  print *, '# ERROR!'
  print *, '# You tried to define the surface', n1, n2, n3, n4
  print *, '# but it doesn''t exists in the block specifications.'
  print *, '# Exiting !'
  stop

1 block=b
  face =fc

  end subroutine
