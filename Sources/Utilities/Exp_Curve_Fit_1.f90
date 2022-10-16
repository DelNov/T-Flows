!==============================================================================!
  program Exp_Curve_Fit
!------------------------------------------------------------------------------!
!   The program fits exponential curve of the form:                            !
!                                                                              !
!   y = a * exp(b*x) + c                                                       !
!                                                                              !
!   through three specified points                                             !
!------------------------------------------------------------------------------!
  implicit none
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: MAX_ITER   = 5000
  integer, parameter :: PLT_POINTS = 1000
  real,    parameter :: TOLERANCE  = 1.0e-8
!-----------------------------------[Locals]-----------------------------------!
  integer :: k
  real    :: x, y, x0, x1, x2, y0, y1, y2
  real    :: dx1, dx2, cx1, cx2, dy1, dy2, dydx1, dydx2
  real    :: y_min, y_max, a_coef, a_coef1, a_coef2, b_coef
  real    :: c_coef, c_coef0, c_coef1, c_coef2, y1c, y2c
  real    :: growth_fac1, growth_fac2, aver_gf, deviation
  real    :: deviation_old, c_coef_fine, ee1, ee2
  real    :: s_x, s_y, s_xy, s_x2, x_avg, y_avg, slope, const
  real    :: y0_ls, y1_ls, y2_ls, deno, sign_coef, delta_c
  real    :: c_coef_old, dev_test, x_scale, delta_c_scale
!==============================================================================!

  !----------------------------------------!
  !                                        !
  !   Hard-code the desired three points   !
  !                                        !
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

  !------------------------------------------!
  !                                          !
  !   Scale x axis in order to avoid large   !
  !    number in expression for deviation    !
  !                                          !
  !------------------------------------------!
  x_scale = max(abs(x0), abs(x1), abs(x2))

  if(x_scale > 10.0) then
    x_scale = 1.0
  else
    x_scale = x_scale * 0.01
  end if

  x0 = x0 / x_scale
  x1 = x1 / x_scale
  x2 = x2 / x_scale

  dx1 = x1 - x0
  dx2 = x2 - x1

  dy1 = y1 - y0
  dy2 = y2 - y1

  cx1 = 0.5*(x0+x1)
  cx2 = 0.5*(x1+x2)

  dydx1 = dy1/dx1
  dydx2 = dy2/dx2

  y_min = min(y0, y1, y2)
  y_max = max(y0, y1, y2)

  !--------------------------------------!
  !                                      !
  !   First estimation of coefficients   !
  !                                      !
  !--------------------------------------!

  ! First estimation of b coef.
  b_coef = log(abs(dydx2/dydx1)) / (cx2 - cx1)

  ! First estimation of a coef.
  a_coef1 = dy1 / (exp(b_coef*x1) - exp(b_coef*x0))
  a_coef2 = dy2 / (exp(b_coef*x2) - exp(b_coef*x1))

  a_coef  = 0.5 * (a_coef1 + a_coef2)

  print *, 'First approx a = ', a_coef
  print *, 'First approx b = ', b_coef

  ! First estimation of c coef.
  c_coef0 = y0 - a_coef * exp(b_coef*x0)
  c_coef1 = y1 - a_coef * exp(b_coef*x1)
  c_coef2 = y2 - a_coef * exp(b_coef*x2)

  c_coef     = (c_coef0 + c_coef1 + c_coef2) / 3.0
  c_coef_old = c_coef

  print *, 'First approx c = ', c_coef

  y1c    = (y1 - c_coef)/(y0 - c_coef)
  y2c    = (y2 - c_coef)/(y1 - c_coef)

  !-------------------------------------------------------!
  !                                                       !
  !   New estimaton of c_coef by forcing that c_coef is   !
  !   not falling into the range defined by y0 and y2.    !
  !                                                       !
  !-------------------------------------------------------!
  c_coef = c_coef/abs(c_coef) * min(abs(y0),abs(y2))

  if(max(abs(dy1/dx1),  &
         abs(dy2/dx2)) / x_scale < 50.0) then
    delta_c_scale = 0.001

  else if(max(abs(dy1/dx1),                         &
              abs(dy2/dx2)) / x_scale > 50.0 .and.  &
          max(abs(dy1/dx1),                         &
              abs(dy2/dx2)) / x_scale < 400) then
    delta_c_scale = 0.00001
  else
    delta_c_scale = 0.000001
  end if

  if(a_coef > 0.0) then
    delta_c = - abs(c_coef) * delta_c_scale
  else
    delta_c =   abs(c_coef) * delta_c_scale
  end if

  c_coef = c_coef + delta_c

  write(*,*) 'New estimation of c_coef ', c_coef

  ! Fine tunning of c coef.
  y1c = (y1 - c_coef) /(y0 - c_coef)
  y2c = (y2 - c_coef) /(y1 - c_coef)

  ee1 = log(y1c)
  ee2 = log(y2c)

  growth_fac1 = exp(ee1/dx1)
  growth_fac2 = exp(ee2/dx2)

  aver_gf = ee1 + ee2

  b_coef  = aver_gf / (x2 - x0)
  write(*,*) 'New estimation of b_coef ', b_coef

  aver_gf = exp(aver_gf/(x2 - x0))

  deviation_old = (abs(aver_gf - growth_fac1)  &
                 + abs(aver_gf - growth_fac2))

  c_coef_old = c_coef

  do k = 1, MAX_ITER
    c_coef_old = c_coef
    c_coef     = c_coef + delta_c

    y1c = ((y1 - c_coef)/(y0 - c_coef))
    y2c = ((y2 - c_coef)/(y1 - c_coef))

    ee1 = log(y1c)
    ee2 = log(y2c)

    growth_fac1 = exp(ee1/dx1)
    growth_fac2 = exp(ee2/dx2)

    aver_gf = ee1 + ee2
    b_coef  = aver_gf / (x2 - x0)

    aver_gf = exp(aver_gf/(x2 - x0))

    deviation = (abs(aver_gf - growth_fac1)   &
               + abs(aver_gf - growth_fac2))

!   print *, 'new dev ', deviation, 'old dev ', deviation_old

    dev_test = deviation_old

    if(deviation < deviation_old) then
      deviation_old = deviation
      c_coef_fine = c_coef
    end if

    if(abs(dev_test - deviation_old) < TOLERANCE) then
      if(k == 1) c_coef_fine = c_coef
      print *, 'Converged in ', k, ' iterations'
      goto 1
    end if
  end do

1 continue

  c_coef = c_coef_fine

  write(*,*) 'c_coef after fine tunning ', c_coef

  !-----------------------------------------------------------------------!
  !                                                                       !
  !   Now we finalize the coefficients by using the least-square method   !
  !                                                                       !
  !-----------------------------------------------------------------------!

  if(c_coef > y_max) then
    y0_ls = log(c_coef - y0)
    y1_ls = log(c_coef - y1)
    y2_ls = log(c_coef - y2)
    sign_coef = -1.0
  else
    y0_ls = log(y0 - c_coef)
    y1_ls = log(y1 - c_coef)
    y2_ls = log(y2 - c_coef)
    sign_coef = 1.0
  end if

  y0_ls = log(abs(y0 - c_coef))
  y1_ls = log(abs(y1 - c_coef))
  y2_ls = log(abs(y2 - c_coef))

  s_x   = x0 + x1 + x2
  s_y   = y0_ls + y1_ls + y2_ls
  s_xy  = x0*y0_ls + x1*y1_ls + x2*y2_ls
  s_x2  = x0**2 + x1**2 + x2**2

  x_avg = s_x / 3.0
  y_avg = s_y / 3.0

  deno = 3.0*s_x2 - s_x**2

  slope = (3.0*s_xy -s_x*s_y) / deno
  const = y_avg - (slope*x_avg)

  a_coef = sign_coef * exp(const)
  b_coef = slope / x_scale
  c_coef = c_coef

  print *, 'a_coef new is ', a_coef
  print *, 'b_coef new is ', b_coef
  print *, 'c_coef new is ', c_coef

  !----------------------------------------!
  !                                        !
  !   Write out the results for checking   !
  !                                        !
  !----------------------------------------!
  x0 = x0 * x_scale
  x1 = x1 * x_scale
  x2 = x2 * x_scale

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
    y = a_coef * exp(b_coef * x) + c_coef
    write(9,*) x, y
  end do
  close(9)
  print *, 'Created curves.dat'

  end program
