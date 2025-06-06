!==============================================================================!
  subroutine get_integer_digits(amg, into, nnum, n_digits, digit)
!------------------------------------------------------------------------------!
!   Decompose non-negative integer into into nnum integers
!
!   This utility is used to parse multi-digit parameters (e.g., ncyc=10110
!   where digits control different solver options). It’s an economic way to
!   encode multiple settings into a single integer.
!
!   Input:  into   - integer (0.le. into .le.999999999)
!           nnum   - integer (1.le. nnum .le.9); number of integers
!                      to be returned on array digit (see below)
!
!   Output: n_digits - integer; number of digits of into
!           digit    - integer-array of length 10:
!                      digit(1)        = first      digit of into,
!                      digit(2)        = second     digit of into,
!                      digit(nnum-1)   = (nnum-1)th digit of into,
!                      digit(nnum)     = rest of into
!                      if nnum > n_digits, the corresponding components
!                      of digit are put to zero.
!
!   Warning: be sure that yout computer can store nnum digits on an
!            integer variable.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(amg_type)      :: amg
  integer, intent(in)  :: into, nnum
  integer, intent(out) :: n_digits, digit(:)
!-----------------------------------[locals]-----------------------------------!
  integer :: i, iq, nrest
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  nrest = into
  n_digits = 1 + int(log10(0.5d0 + dble(into)))

  if(nnum .lt. n_digits) then
    iq = 10**(n_digits - nnum + 1)
    digit(nnum) = nrest-(nrest/iq)*iq
    nrest = nrest/iq
    do i = nnum-1, 1, -1
      digit(i) = nrest-(nrest/10)*10
      nrest = nrest/10
    end do
    return
  end if

  do i = n_digits, 1, -1
    digit(i) = nrest-(nrest/10)*10
    nrest = nrest/10
  end do

  do i = n_digits+1, nnum
    digit(i) = 0
  end do

  end subroutine
