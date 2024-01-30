!==============================================================================!
  program Primes
!------------------------------------------------------------------------------!
!   This program finds all prime numbers in the given last                     !
!                                                                              !
!   Compile with: gfortran -o Primes Primes.f90                                !
!------------------------------------------------------------------------------!
  implicit  none
!------------------------------------------------------------------------------!
  integer :: first, last, number, divisor, count
!==============================================================================!

  ! Keep trying to read a good input
  do
    print *,  'At which number you want to begin?';  read(*,*) first
    print *,  'At which number you want to end?';    read(*,*) last

    if(first .gt. 1 .and. last .gt. first) exit

    if(first .le. 1) then
      print *, 'The range should start at number 2 or bigger number!'
    end if

    if(last .le. first) then
      print *, 'Come on! The end of range should be larger than the beginning!'
    end if

    print *, 'Please try again!'
    print *, ''
  end do

  ! Input is correct; start counting
  count = 0

  ! If first is even, increase it by one, to make it odd
  if(first .ne. 2 .and. mod(first,2) .eq. 0) first = first + 1

  ! Number 2 is kind of a special case, deal with it
  if(first .eq. 2) then
    count = count + 1                    ! yes, two is a prime number
    print *,  'Prime number #', count, ': ', 2
    first = 3
  end if

  do number = first, last, 2             ! try all odd numbers

    divisor = 3                          ! divisor starts with 3
    do
      if (divisor*divisor > number .or. mod(number,divisor) == 0)  exit
      divisor = divisor + 2              ! if does not evenly divide, next odd
    end do

    if (divisor*divisor > number) then   ! are all divisor exhausted?
      count = count + 1                  ! yes, this number is a prime
      print *,  'Prime number #', count, ': ', number
    end if
  end do

  print *, ''
  print *, 'There are ', count, ' primes in the last of 2 and ', last

  end program
