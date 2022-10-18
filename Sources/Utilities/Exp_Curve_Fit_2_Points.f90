  include 'Const_Mod.f90'
  include 'Math_Mod.f90'

!==============================================================================!
  program Exp_Curve_Fit
!------------------------------------------------------------------------------!
!   The program fits exponential curve of the form:                            !
!                                                                              !
!   y = a * exp(b*x) + c                                                       !
!                                                                              !
!------------------------------------------------------------------------------!
  use Math_Mod
!------------------------------------------------------------------------------!
  implicit none
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: PLT_POINTS = 1000
!-----------------------------------[Locals]-----------------------------------!
  real    :: x, y, x0, y0, x1, y1, x2, y2, a, b, c, dy_dx_0
  integer :: k
!==============================================================================!

  !----------------------------------!
  !   Hard-code the desired values   !
  !----------------------------------!
  x0 = 0.0
  y0 = 0.5

  x1 = 0.005
  y1 = -5.3244798523939757E-004

  x2 = 0.015
  y2 = -8.8180312077764597E-002

  dy_dx_0 = -222.366547   ! value obtained from other program

  print *, '--------------'
  print *, 'Initial values'
  print *, '--------------'
  print *, 'dy/dx at 0 : ', dy_dx_0
  print *, 'x0, y0     : ', x0, 'unknown'
  print *, 'x1, y1     : ', x1, y1
  print *, 'x2, y2     : ', x2, y2

  !-------------------------------------------------------------!
  !   Fit exponential function with derivative and two points   !
  !-------------------------------------------------------------!
  call Math % Fit_Exp_Derivative_And_Two_Points(dy_dx_0,  &
                                                x0, y0,   &
                                                x1, y1,   &
                                                x2, y2,   &
                                                a, b, c)

  print *, '---------------'
  print *, 'Computed values'
  print *, '---------------'
  print *, 'a_coef        : ', a
  print *, 'b_coef        : ', b
  print *, 'c_coef        : ', c
  print *, 'dy/dx at 0 is ', dy_dx_0
  print *, 'x0, y0     : ', x0, y0
  print *, 'x1, y1     : ', x1, y1
  print *, 'x2, y2     : ', x2, y2

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
