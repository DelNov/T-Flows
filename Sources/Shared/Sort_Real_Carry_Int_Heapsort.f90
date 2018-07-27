!==============================================================================!
  subroutine Sort_Real_Carry_Int_Heapsort(a, b, n)
!------------------------------------------------------------------------------!
! taken from                                                                   !
! https://rosettacode.org/wiki/Sorting_algorithms/Heapsort#Fortran             !
! https://en.wikipedia.org/wiki/Heapsort                                       !
! modified to carry int array also                                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer                 :: n
  real,    intent(in out) :: a(0: n-1)
  integer, intent(in out) :: b(0: n-1)
!-----------------------------------[Locals]-----------------------------------!
  integer                 :: start, bottom
!==============================================================================!
 
  do start = (n - 2) / 2, 0, -1
    call siftdown(a, start, n)
  end do

  do bottom = n - 1, 1, -1
    call Swap_Real(a(0), a(bottom))
    call Swap_Int (b(0), b(bottom))
    call siftdown(a, 0, bottom)
  end do
!------------------------------------------------------------------------------!
contains
!------------------------------------------------------------------------------!
subroutine siftdown(a, start, bottom)
!------------------------------------------------------------------------------!
  real, intent(in out) :: a(0:)
  integer, intent(in)  :: start, bottom
  integer              :: child, root
!------------------------------------------------------------------------------!
 
  root = start
  do while(root*2 + 1 < bottom)
    child = root * 2 + 1
 
    if (child + 1 < bottom) then
      if (a(child) < a(child+1)) child = child + 1
    end if
 
    if (a(root) < a(child)) then
      call Swap_Real(a(child), a(root))
      call Swap_Int (b(child), b(root))
      root = child
    else
      return
    end if  
  end do      
 
end subroutine siftdown
!------------------------------------------------------------------------------!
 
end subroutine