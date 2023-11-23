!==============================================================================!
  subroutine Fit_Exp_Three_Points(Math,     &
                                  dy_dx_0,  &  ! unknown
                                  x0, y0,   &  ! x0=0, y0 is prescribed
                                  x1, y1,   &  ! prescribed
                                  x2, y2)      ! prescribed
!------------------------------------------------------------------------------!
!>  This subroutine fits exponential curve of the form:
!>  y = a * exp(b*x) + c
!>  through three specified points.
!------------------------------------------------------------------------------!
!   This subroutine fits exponential curve of the form:                        !
!                                                                              !
!   y = a * exp(b*x) + c                                                   (1) !
!                                                                              !
!   through three specified points.                                            !
!   It is assumed, however, that coordinate system starts at zero.  This is    !!
!   not too much of a nuissance, since "wall coordinate" starts at zero.       !
!                                                                              !
!   The solutuion works as follows. Since x0 is zero, at x0:                   !
!                                                                              !
!   y0 = a + c,  or:                                                       (2) !
!   c  = y0 - a                                                            (3) !
!                                                                              !
!   So, in general, it is valid that:                                          !
!                                                                              !
!   y = a * exp(b*x) + y0 - a                                              (4) !
!   y = y0 + a * (exp(b*x) - 1)                                            (5) !
!                                                                              !
!   If we define equation (5) for points 1 and 2, we get:                      !
!                                                                              !
!   y1 = y0 + a * (exp(b*x1) - 1)                                          (6) !
!   y2 = y0 + a * (exp(b*x2) - 1)                                          (7) !
!                                                                              !
!   or:                                                                        !
!                                                                              !
!   y1 - y0 = a * (exp(b*x1) - 1)                                          (8) !
!   y2 - y0 = a * (exp(b*x2) - 1)                                          (9) !
!                                                                              !
!   Dividing equations (8) and (9) yields:                                     !
!                                                                              !
!   y1 - y0   exp(b*x1) - 1                                                    !
!   ------- = -------------                                               (10) !
!   y2 - y0   exp(b*x2) - 1                                                    !
!                                                                              !
!        exp(b*x1) - 1              y1 - y0                                    !
!   d =  -------------  where:  d = -------                               (11) !
!        exp(b*x2) - 1              y2 - y0                                    !
!                                                                              !
!   It is now about the time we come to the point.  This subroutine solves     !
!   equation (11) in a Newton-Raphson like fashion to obtain coefficient "b".  !
!   Once "b" is known, coefficient "a" is obatined from equation (7) and       !
!   coefficient "c" from equation (3).                                         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Math_Type)  :: Math
  real, intent(out) :: dy_dx_0  !! unknown derivative dy/dx at 0
  real, intent(in)  :: x0, y0   !! x0 is zero, y0 prescribed
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
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(x0)  ! x0 is kept for symmetry with sister procedure
!==============================================================================!

  !-----------------------------!
  !   Calculate coefficient d   !
  !-----------------------------!
  d = (y1 - y0) / Math % Signed_Lower_Limit((y2 - y0), PICO)

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
      e_array(k) = (exp(b_array(k)*x1) - 1.0)                                 &
                 / Math % Signed_Lower_Limit(exp(b_array(k)*x2) - 1.0, PICO)  &
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
  !   At this point, coefficient b is known and we can compute "a" from (9)   !
  !---------------------------------------------------------------------------!
  a = (y2 - y0) / Math % Signed_Lower_Limit(exp(b*x2) - 1.0, PICO)

  !------------------------------------------------!
  !   Coefficient "c" can be computed from (3)     !
  !------------------------------------------------!
  c = y0 - a

  !---------------------------------------!
  !   Calculate remaining unkown values   !
  !---------------------------------------!
  dy_dx_0 = a * b

  end subroutine
