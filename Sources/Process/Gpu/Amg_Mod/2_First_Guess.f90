!==============================================================================!
  subroutine First_Guess(Amg, ifirst, u)
!------------------------------------------------------------------------------!
!   Puts a first approximation to finest grid
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(Amg_Type) :: Amg
  integer         :: ifirst
  real            :: u(:)
!-----------------------------------[locals]-----------------------------------!
  real    :: s, sd
  integer :: digit(AMG_MAX_LEVELS)
  integer :: i, ifrst, n_digits
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  if(ifirst .ne. 0) then
    call Amg % Get_Integer_Digits(ifirst, 3, n_digits, digit)
    ifrst = digit(2)
  else
    ifrst = 3
  endif

  if(ifrst.eq.1) then
    do i = Amg % imin(1), Amg % imax(1)
      u(i) = 0.0
    end do
    return

  else if(ifrst.eq.2) then
    do i = Amg % imin(1), Amg % imax(1)
      u(i) = 1.0
    end do
    if(Amg % irow0 .eq. AMG_SINGULAR_MATRIX) u(Amg % imax(1)) = 0.0
    return

  else if(ifrst.eq.3) then
    sd = 0.72815
    if(digit(3)*ifirst .ne. 0) then
      s = real(digit(3))
      do i = 1, 10
        s = s * 0.1
        if(s .lt. 1.0) exit
      end do
    else
      s = sd
    end if
    do i = Amg % imin(1), Amg % imax(1)
      u(i) = Amg % Random_0_To_0p1(s)
    end do
    if(Amg % irow0 .eq. AMG_SINGULAR_MATRIX) u(Amg % imax(1)) = 0.0
    return
  endif

  end subroutine
