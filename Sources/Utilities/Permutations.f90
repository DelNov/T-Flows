
subroutine Swap_Int(a, b)

  integer :: a, b
  integer :: t

  t = a
  a = b
  b = t

end subroutine

recursive subroutine Generate(k, A)

  integer :: k
  integer :: A(:)

  if(k .eq. 1) then
    print *, A
  else
    ! Generate permutations with kth unaltered
    ! Initially k == length(A)
    call Generate(k - 1, A)

    ! Generate permutations for kth swapped with each k-1 initial
    do i = 1, k - 1
      ! Swap choice dependent on parity of k (even or odd)
      if(mod(k,2) .eq. 0) then
        call Swap_Int(A(i), A(k)) !  zero-indexed, the kth is at k-1
      else
        call Swap_Int(A(1), A(k))
     end if
     call Generate(k - 1, A)

    end do
  end if

end subroutine

program Caller

  interface
    recursive subroutine Generate(k, A)
    integer :: k
    integer :: A(:)
    end subroutine
  end interface


  integer :: k

  integer :: A(3) = (/1, 2, 3/)

  call Generate(3, A)

end program

