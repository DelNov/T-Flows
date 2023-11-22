!==============================================================================!
  integer function Is_Line_In_Block(Dom, n1, n2, b)
!------------------------------------------------------------------------------!
!   Checks if the line defined n1 and n2 is inside the block b.                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Domain_Type)  :: Dom
  integer, intent(in) :: b, n1, n2
!-----------------------------------[Locals]-----------------------------------!
  integer :: l1, l2
!==============================================================================!

  do l1=1,8
    do l2=1,8
      if( (Dom % blocks(b) % corners(l1) .eq. n1) .and.   &
          (Dom % blocks(b) % corners(l2) .eq. n2) ) then
           goto 1
      end if
    end do
  end do

  Is_Line_In_Block = 0
  return

1 Is_Line_In_Block = b
  return

  end function
