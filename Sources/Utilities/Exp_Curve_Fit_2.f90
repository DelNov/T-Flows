!==============================================================================!
  subroutine Set_Range(n, minv, maxv, array)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer, intent(in)  :: n
  real   , intent(in)  :: minv, maxv
  real,    intent(out) :: array(n)
!-----------------------------------[Locals]-----------------------------------!
  integer :: i
  real    :: delta
!==============================================================================!

  array(1) = minv
  array(n) = maxv
  delta = (maxv-minv) / real(n-1.0)
  do i = 1, n
    array(i) = minv + delta * (i-1)
  end do

  end subroutine

!==============================================================================!
  program Exp_Curve_Fit
!------------------------------------------------------------------------------!
!   The program fits exponential curve of the form:                            !
!                                                                              !
!   y = a * exp(b*x) + c                                                       !
!                                                                              !
!   through three specified points.  It is assumed, however, that coordinate   !
!   system starts at zero.  Like, the wall starts at zero.                     !
!------------------------------------------------------------------------------!
  implicit none
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: MAX_ITER   =   64
  integer, parameter :: PLT_POINTS = 1024
  integer, parameter :: N_SAMPLES  =    8
  real,    parameter :: RANGE_MIN  =   -1.0e+5
  real,    parameter :: RANGE_MAX  =   +1.0e+5
  real,    parameter :: TOLERANCE  =    1.0e-8
  logical, parameter :: DEBUG      = .true.
!-----------------------------------[Locals]-----------------------------------!
  integer :: i, j, k
  real    :: x, y, x0, y0, x1, y1, x2, y2
  real    :: a, b, c, d, b_array(N_SAMPLES), e_array(N_SAMPLES)
!==============================================================================!

  !----------------------------------------!
  !   Hard-code the desired three points   !
  !----------------------------------------!
  x0 = 0.0
  y0 = 0.5

  x1 = 0.005
  y1 = -5.3244798523939757E-004

  x2 = 0.015
  y2 = -8.8180312077764597E-002

  print *, 'x0, y0, ', x0, y0
  print *, 'x1, y1, ', x1, y1
  print *, 'x2, y2, ', x2, y2

  !------------------------------------------------!
  !   Since x0 is zero, at x0:                     !
  !                                                !
  !   y0 = a + c  or                               !
  !   c  = y0 - a                                  !
  !   y = a * exp(b*x) + y0 - a                    !
  !   y = y0 + a * (exp(b*x) - 1)                  !
  !                                                !
  !   Dividing equations at x1 and x2:             !
  !                                                !
  !   y1 = y0 + a * (exp(b*x1) - 1)                !
  !   y2 = y0 + a * (exp(b*x2) - 1)                !
  !                                                !
  !   y1 - y0 = a * (exp(b*x1) - 1)                !
  !   y2 - y0 = a * (exp(b*x2) - 1)                !
  !                                                !
  !   we get:                                      !
  !                                                !
  !   y1 - y0   exp(b*x1) - 1                      !
  !   ------- = -------------                      !
  !   y2 - y0   exp(b*x2) - 1                      !
  !                                                !
  !   or:                                          !
  !        exp(b*x1) - 1              y1 - y0      !
  !   d =  -------------  where:  d = -------      !
  !        exp(b*x2) - 1              y2 - y0      !
  !                                                !
  !   The last expression will be used in Newton   !
  !   Raphson like iterations (within samples) to  !
  !   minimize the error:                          !
  !                                                !
  !        exp(b*x1) - 1                           !
  !   e =  -------------  - d                      !
  !        exp(b*x2) - 1                           !
  !------------------------------------------------!

  ! Calculate coefficient d
  d = (y1 - y0) / (y2 - y0)

  if(DEBUG) print *, 'd = ', d

  ! Set initial range for b
  call Set_Range(N_SAMPLES, RANGE_MIN, RANGE_MAX, b_array)

  !--------------------------------------------------!
  !                                                  !
  !   Browse throut Newton-Raphson like iterations   !
  !                                                  !
  !--------------------------------------------------!
  do i = 1, MAX_ITER

    ! Compute errors for the entire range of samples
    do k = 1, N_SAMPLES
      e_array(k) = (exp(b_array(k)*x1) - 1.0)   &
                 / (exp(b_array(k)*x2) - 1.0)   &
                 - d
    end do

    ! Find the two samples which change the sign in error and set
    ! them as new bounding values for array of coefficitens b
    do j = 2, N_SAMPLES
      if(e_array(j-1) * e_array(j) < 0.0) then
        call Set_Range(N_SAMPLES, b_array(j-1), b_array(j), b_array)
        goto 1
      end if
    end do
1   continue

    ! Check if desired tolerance has been reached
    if(abs(b_array(1) - b_array(N_SAMPLES)) < TOLERANCE) goto 2

  end do

2 continue
  if(DEBUG) print *, 'Converged in ', i, 'iterations'
  b = 0.5 * (b_array(1) + b_array(N_SAMPLES))

  !-------------------------------------------!
  !   At this point, coefficient b is known   !
  !   and we can compute a from:              !
  !                                           !
  !   y2 - y0 = a * (exp(b*x2) - 1)           !
  !                                           !
  !   or:                                     !
  !           y2 - y0                         !
  !   a = ---------------                     !
  !       (exp(b*x2) - 1)                     !
  !                                           !
  !-------------------------------------------!
  a = (y2 - y0) / (exp(b*x2) - 1.0)

  !------------------------------------------------!
  !   We established before that:                  !
  !                                                !
  !   y0 = a + c  or:                              !
  !                                                !
  !   c = y0 - a                                   !
  !------------------------------------------------!
  c = y0 - a

  !----------------------------------------------!
  !                                              !
  !   At this point, coefficients are computed   !
  !                                              !
  !----------------------------------------------!
  print *, 'a_coef new is ', a
  print *, 'b_coef new is ', b
  print *, 'c_coef new is ', c

  !-----------------------------------!
  !   Write out results for xmgrace   !
  !-----------------------------------!
  open(9, file='points.dat')
  write(9,*) x0, y0
  write(9,*) x1, y1
  write(9,*) x2, y2
  close(9)
  print *, 'Created points.dat'

  !-----------------------------------!
  !   Write out results for xmgrace   !
  !-----------------------------------!
  open(9, file='curve.dat')
  do k = 0, PLT_POINTS
    x = x0 + (x2-x0)/real(PLT_POINTS) * k
    y = a * exp(b * x) + c
    write(9,*) x, y
  end do
  close(9)
  print *, 'Created curves.dat'

  end program
