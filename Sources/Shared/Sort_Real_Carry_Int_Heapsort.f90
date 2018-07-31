!==============================================================================!
  subroutine Sort_Real_Carry_Int_Heapsort(a, b, n)
!------------------------------------------------------------------------------!
!   Taken from:                                                                !
!   https://rosettacode.org/wiki/Sorting_algorithms/Heapsort#Fortran           !
!                                                                              !
!   Algorith is described in:                                                  !
!   https://en.wikipedia.org/wiki/Heapsort                                     !
!   modified to carry int array b, to accept arrays starting from 1:           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer                 :: n
  real,    intent(in out) :: a(1: n)
  integer, intent(in out) :: b(1: n)
!-----------------------------------[Locals]-----------------------------------!
  integer                 :: start, bottom
!==============================================================================!

  do start = (n - 2) / 2 + 1, 1, -1
    call Sift_Down(a, start, n + 1)
  end do

  do bottom = n, 2, -1
    call Swap_Real(a(1), a(bottom))
    call Swap_Int (b(1), b(bottom))
    call Sift_Down(a, 1, bottom)
  end do

  contains

  !============================================================================!
    subroutine Sift_Down(a, start, bottom)
  !----------------------------------------------------------------------------!
  !--------------------------------[Arguments]---------------------------------!
    real, intent(in out) :: a(1:)
    integer, intent(in)  :: start, bottom
  !----------------------------------[Locals]----------------------------------!
    integer              :: child, root
  !============================================================================!

    root = start

    do while( (root-1) * 2 + 2 < bottom)
      child = (root-1) * 2 + 2

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

  end subroutine

end subroutine