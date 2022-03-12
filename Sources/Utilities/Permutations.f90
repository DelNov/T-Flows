!==============================================================================!
  subroutine Swap_Int(a, b)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: a, b
!-----------------------------------[Locals]-----------------------------------!
  integer :: t
!==============================================================================!

  t = a
  a = b
  b = t

  end subroutine

!==============================================================================!
  recursive subroutine Generate(k, a)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: k
  integer :: a(:)
!-----------------------------------[Locals]-----------------------------------!
  integer :: i
!==============================================================================!

  if(k .eq. 1) then
    print *, a
  else
    ! Generate permutations with kth unaltered
    ! Initially k == length(a)
    call Generate(k - 1, a)

    ! Generate permutations for kth swapped with each k-1 initial
    do i = 1, k - 1
      ! Swap choice dependent on parity of k (even or odd)
      if(mod(k,2) .eq. 0) then
        call Swap_Int(a(i), a(k)) !  zero-indexed, the kth is at k-1
      else
        call Swap_Int(a(1), a(k))
     end if
     call Generate(k - 1, a)

    end do
  end if

  end subroutine

!==============================================================================!
  program Caller
!----------------------------------[Calling]-----------------------------------!
  interface
    recursive subroutine Generate(k, a)
    integer :: k
    integer :: a(:)
    end subroutine
  end interface
!-----------------------------------[Locals]-----------------------------------!
  integer              :: i, k
  integer, allocatable :: a(:)
  character(80)        :: arg
!==============================================================================!

  ! Proper invocation
  if(command_argument_count() .eq. 1) then

    ! Get k
    call get_command_argument(1, arg);  read(arg, *) k

    allocate(a(k))
    do i = 1, k
      a(i) = i
    end do

    call Generate(k, a)

  ! Wrong invocation
  else
    print *, '# You failed to invoke the program properly.'
    print *, '# Correct invocation:'
    print *, './Permutations  number'
    print *, '# where "number" is the range for permutations'
    stop
  end if

  end program

