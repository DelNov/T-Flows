!==============================================================================!
  subroutine first_guess(amg, ifirst, u)
!------------------------------------------------------------------------------!
!   Puts a first approximation to finest grid
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(amg_type)  :: amg
  integer          :: ifirst
  double precision :: u(:)
!-----------------------------------[locals]-----------------------------------!
  double precision :: s, sd
  integer          :: digit(AMG_MAX_LEVELS)
  integer          :: i, ifrst, n_digits
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  if(ifirst .ne. 0) then
    call amg % get_integer_digits(ifirst, 3, n_digits, digit)
    ifrst = digit(2)
  else
    ifrst = 3
  endif

  if(ifrst.eq.1) then
    do i = amg % imin(1), amg % imax(1)
      u(i) = 0.0d0
    end do
    return

  else if(ifrst.eq.2) then
    do i = amg % imin(1), amg % imax(1)
      u(i) = 1.0d0
    end do
    if(amg % irow0 .eq. AMG_SINGULAR_MATRIX) u(amg % imax(1)) = 0.0d0
    return

  else if(ifrst.eq.3) then
    sd = 0.72815d0
    if(digit(3)*ifirst .ne. 0) then
      s = dble(digit(3))
      do i = 1, 10
        s = s*0.1d0
        if(s .lt. 1.0d0) exit
      end do
    else
      s = sd
    end if
    do i = amg % imin(1), amg % imax(1)
      u(i) = amg % random_0_to_0p1(s)
    end do
    if(amg % irow0 .eq. AMG_SINGULAR_MATRIX) u(amg % imax(1)) = 0.0d0
    return
  endif

  end subroutine
