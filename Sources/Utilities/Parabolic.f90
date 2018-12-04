!==============================================================================!
  program Parabolic
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  real    :: x_s =  -0.005
  real    :: x_e =   0.005
  real    :: bulk = -1.0
  integer :: i
  integer :: n = 41

  real    :: x_l, x_r, x_c, y_l, y_r, y_c, dx, lx, integral, area
!==============================================================================!

  print *, '#=================='
  print *, '# Number of points '
  print *, '#=================='
  print *, n

  print *, '#================================='
  print *, '#   Coordinate      Velocity'
  print *, '#================================='

  ! distance between two nodes
  lx = (x_e - x_s)
  dx = lx / (n-1)

  ! Print profile and compute integral
  area     = 0.0
  integral = 0.0
  do i = 1, n
    x_l  = (x_s + dx * (i-1))  ! x left
    x_r  = (x_s + dx *  i   )  ! x rigth
    y_l = (1.0 - ((x_l/lx) ** 2 / 0.25)) * bulk *  1.5
    y_r = (1.0 - ((x_r/lx) ** 2 / 0.25)) * bulk *  1.5

    print '(2es16.5)', x_l, y_l

    if(i < n) then
      y_c      = (y_l + y_r) * 0.5
      area     = area     + dx
      integral = integral + y_c * dx
    end if
  end do

  print *, '# Bulk velocity is: ', integral / area

  ! Compute integral

  end program
