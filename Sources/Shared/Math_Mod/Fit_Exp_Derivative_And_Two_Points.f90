!==============================================================================!
  subroutine Fit_Exp_Derivative_And_Two_Points(Math,     &
                                               dy_dx_0,  &  ! prescribed
                                               x0, y0,   &  ! x0=0, y0 unknown
                                               x1, y1,   &  ! prescribed
                                               x2, y2)      ! prescribed
!------------------------------------------------------------------------------!
!>  This subroutine fits exponential curve of the form:
!>  y = a * exp(b*x) + c
!>  through three specified points, where first point starts at zero.
!>  In the first point, a derivative (in x) is given, and the values are
!>  given at the remining points.
!------------------------------------------------------------------------------!
!   This subroutine fits exponential curve of the form:                        !
!                                                                              !
!   y = a * exp(b*x) + c                                                   (1) !
!                                                                              !
!   through three specified points, where first point starts at zero.          !
!   In the first point, a derivative (in x) is given, and the values are       !!
!   given at the remining points.                                              !
!                                                                              !
!   The solutuion works as follows. Derivative of (1) with respect to x is:    !
!                                                                              !
!   dy/dx = a * b * exp(b*x)                                               (2) !
!                                                                              !
!   At x = 0:                                                                  !
!                                                                              !
!   (dy/dx)_0 = a * b                                                      (3) !
!                                                                              !
!   Moreover, for points 1 and 2 (where x is non-zero), we have:               !
!                                                                              !
!   y1 = a * exp(b*x1) + c                                                 (4) !
!                                                                              !
!   y2 = a * exp(b*x2) + c                                                 (5) !
!                                                                              !
!   Substracting (5) from (4):                                                 !
!                                                                              !
!   y2 - y1 = a * (exp(b*x2) - exp(b*x1))                                  (6) !
!                                                                              !
!   and inserting (3), we get:                                                 !
!                                                                              !
!                       exp(b*x2) - exp(b*x1)                                  !
!   y2 - y1 = (dy/dx)_0 ---------------------                              (7) !
!                                 b                                            !
!                                                                              !
!   or:                                                                        !
!                                                                              !
!    y2 - y1    exp(b*x2) - exp(b*x1)                                          !
!   --------- = ---------------------                                      (8) !
!   (dy/dx)_0             b                                                    !
!                                                                              !
!   If we define d as:                                                         !
!                                                                              !
!        y2 - y1                                                               !
!   d = ---------                                                          (9) !
!       (dy/dx)_0                                                              !
!                                                                              !
!   equation (8) assumes the form:                                             !
!                                                                              !
!       exp(b*x2) - exp(b*x1)                                                  !
!   d = ---------------------                                             (10) !
!                 b                                                            !
!                                                                              !
!   It is now about the time we come to the point.  This subroutine solves     !
!   equation (9) in a Newton-Raphson like fashion to obtain coefficient "b".   !
!   Once "b" is known, coefficient "a" is obatined from equation (3) and       !
!   coefficient "c" from equation (5).                                         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Math_Type)  :: Math
  real, intent(in)  :: dy_dx_0  !! prescribed dy/dx at 0
  real, intent(out) :: x0, y0   !! x0 is zero, y0 is unknown
  real, intent(in)  :: x1, y1   !! prescribed point 1
  real, intent(in)  :: x2, y2   !! prescribed point 2
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: MAX_ITER  = 64
  integer, parameter :: N_SAMPLES =  8
  real,    parameter :: REL_TOLER =  1.0e-3
  logical, parameter :: DEBUG     = .false.
!-----------------------------------[Locals]-----------------------------------!
  integer :: i, j, k
  real    :: a, b, c, d, b_array(N_SAMPLES), e_array(N_SAMPLES)
!==============================================================================!

  !-----------------------------!
  !   Calculate coefficient d   !
  !-----------------------------!
  d = (y2 - y1) / Math % Signed_Lower_Limit(dy_dx_0, PICO)

  ! Set initial range for b
  call Math % Set_Array_Range(N_SAMPLES, -100.0/abs(x2),  &
                                         +100.0/abs(x2),  &
                                          b_array)

  !---------------------------------------------------------------!
  !   Browse through Newton-Raphson like iterations to find "b"   !
  !---------------------------------------------------------------!
  do i = 1, MAX_ITER

    if(DEBUG) print *, '---------------------------------'
    if(DEBUG) print *, 'Iteration ', i
    if(DEBUG) print *, '---------------------------------'

    ! Compute errors for the entire range of samples
    do k = 1, N_SAMPLES
      e_array(k) = (exp(b_array(k)*x2) - exp(b_array(k)*x1))    &
                 / Math % Signed_Lower_Limit(b_array(k), PICO)  &
                 - d
      if(DEBUG) print *, 'e_array(k) = ', e_array(k)
    end do

    ! Find the two samples which change the sign in error and set
    ! them as new bounding values for array of coefficitens b.
    !
    ! In other words: "zoom in" the sector of b_array where error
    ! changed the sign for the next iteration.
    do j = 2, N_SAMPLES
      if(e_array(j-1) * e_array(j) < 0.0) then
        call Math % Set_Array_Range(N_SAMPLES, b_array(j-1),  &
                                               b_array(j),    &
                                               b_array)
        goto 1
      end if
    end do
1   continue

    ! Check if desired relative tolerance has been reached
    if(  abs(b_array(1) - b_array(N_SAMPLES))  &
       < abs(b_array(1) + b_array(N_SAMPLES)) * REL_TOLER) goto 2

  end do

  !---------------------------------------------------------!
  !   This is the point at which coefficent "b" converged   !
  !---------------------------------------------------------!
2 continue
  b = 0.5 * (b_array(1) + b_array(N_SAMPLES))

  !---------------------------------------------------------------------------!
  !   At this point, coefficient b is known and we can compute "a" from (3)   !
  !---------------------------------------------------------------------------!
  a = dy_dx_0 / Math % Signed_Lower_Limit(b, PICO)

  !------------------------------------------------!
  !   Coefficient "c" can be computed from (5)     !
  !------------------------------------------------!
  c = y2 - a * exp(b*x2)

  !---------------------------------------!
  !   Calculate remaining unkown values   !
  !---------------------------------------!
  x0 = 0
  y0 = a + c

  end subroutine
