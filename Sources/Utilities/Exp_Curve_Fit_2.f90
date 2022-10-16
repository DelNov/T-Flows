!==============================================================================!
  program Exp_Curve_Fit
!------------------------------------------------------------------------------!
!   The program fits exponential curve of the form:                            !
!                                                                              !
!   y = a * exp(b*x) + c                                                       !
!                                                                              !
!   Compile with:                                                              !
!   gfortran -c ../Shared/Math_Mod.f90                                         !
!   gfortran -o exp_fit Math_Mod.o Exp_Curve_Fit_2.f90                         !
!------------------------------------------------------------------------------!
  use Math_Mod
!------------------------------------------------------------------------------!
  implicit none
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: PLT_POINTS = 1000
!-----------------------------------[Locals]-----------------------------------!
  real    :: x, y, x0, y0, x1, y1, x2, y2, a, b, c
  integer :: k
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

  !---------------------------------------------------!
  !   Fit Exponential Function Through Three Points   !
  !---------------------------------------------------!
  call Math % Fit_Exponential_Function_Through_Three_Points(x0, y0,  &
                                                            x1, y1,  &
                                                            x2, y2,  &
                                                            a, b, c)

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
